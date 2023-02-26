setwd("~/Documents/FW")
library(tidyverse)

## import data
data <- read_tsv("Meta_Counts_T.txt", comment="#")
data <- data[!duplicated(data$Compound), ]
SampleID <- c( "H1", "H2", "H3", "H4", "H5", "C1", "C2", "C3","C4", "C5")
Group <-c(rep("H", 5),rep("C", 5))


## reformat data
countdata <- data %>%
  as.data.frame() %>% 
  column_to_rownames("Compound") %>% # turn the geneid column into rownames
  as.matrix()

logdata <- log(countdata, 2)%>% as.data.frame()


paretoscale <- function(z) {
  rowmean <- apply(z, 1, mean)
  rowsd <- apply(z, 1, sd)
  rowsqrtsd <- sqrt(rowsd)
  rv <- sweep(z, 1, rowmean, "-")
  rv <- sweep(rv, 1, rowsqrtsd, "/")
  return(rv)
}
  
pareto.logdata <- paretoscale(logdata) %>% as.data.frame()


calc_ttest <- function(df, gr1, gr2, maxAdjP, minFC) {
  index_g1 <- which(Group %in% c(gr1))
  index_g2 <- which(Group %in% c(gr2))
    #Log2 fold change group2 - group1
  df$Log2FC <- apply(
    df, 1, function(x) {
      mean( x[index_g1] ) - mean( x[index_g2] )
    }
  )
  
  #T-test with equal variance
  df$T_Pval <- apply(
    df, 1, function(x) {
      res <- t.test(
        x[index_g1], x[index_g2],
        alternative = "two.sided", var.equal = TRUE
      )
      mean(x[index_g1]) - mean(x[index_g2] )
      res$p.value
    }
  )
  #Benjamini-Hochberg correction for multiple testing
  df$adjPval <- p.adjust(df$T_Pval, method = "BH")
  #Add the categorical column for easier visualization
  df$Diff_Abund <- apply(
    df, 1, function(x) {
      if (x[["adjPval"]] <= maxAdjP & x[["Log2FC"]] >= minFC) {
        return( paste("Up in", gr1) )
      } else if (x[["adjPval"]] <= maxAdjP & x[["Log2FC"]] <= -1*minFC) {
        return( paste("Up in", gr2) )
      } else {
        return('Non-significant')
      }
    }
  )
  df
}
maxAdjP <- 0.2
minLog2FC <- 0
gr1 <- "H"
gr2 <- "C"
dfTtest <- calc_ttest(logdata, gr1, gr2, maxAdjP, minLog2FC)

write.table(dfTtest, file="metabolites_T.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)


sampleinfo_input <- data.frame(SampleID, Group)%>% 
  column_to_rownames("SampleID")

MyResult.splsda <- splsda(t(logdata), sampleinfo$Group, keepX = c(50,50)) 
plotIndiv(MyResult.splsda)
plotVar(MyResult.splsda)
