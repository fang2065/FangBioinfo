library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(limma)
library(qvalue)
library(ggfortify)
library(ggrepel)
library(UniProt.ws)
library(DOSE)
library(org.Mm.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(ggridges)
library(GseaVis) #https://github.com/junjunlab/GseaVis
library(mixOmics) 

data <- read_tsv("Normalised_peptides_counts.txt", comment="#")
data<-data[!is.na(data$H1) & !is.na(data$C1),] #remove NA value 

countdata <- data %>%
  as.data.frame() %>% 
  column_to_rownames("Accession") # turn the peptide column into rownames



#box plot 
par(mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
boxplot(countdata[, 1:6],  ylim = c(0, 2.5), main="Boxplot normalized Intensities")

sumstats <- function(count) {
    Mean <- apply(count, 1, mean)
    Median <- apply(count, 1, median)
    SD <- apply(count, 1, sd)
    SE <- apply(count, 1, function(x) sd(x)/sqrt(length(x)))
    CV <- apply(count, 1, function(x) sd(x)/mean(x))
    result <- data.frame(Mean, Median, SD, SE, CV)
    return(result)
}

sumstats_input <- sumstats(countdata)

#log2 transformation if necessary 
#countdata <- log2(countdata)

#PCA
MyResult.pca <- PCA(countdata)
plotIndiv(MyResult.pca)
#or sPLS-DA
SampleID <- c("H1", "H2", "H3", "C1", "C2", "C3")
Group <-c(rep("C", 3),rep("H", 3))
sampleinfo <- data.frame(SampleID, Group)%>% 
  column_to_rownames("SampleID") # generate sampleinfo file 
MyResult.splsda <- splsda(countdata, sampleinfo, keepX = c(50,50)) 
plotIndiv(MyResult.splsda)

#two group t test
calc_ttest <- function(df, groupping, gr1, gr2, maxAdjP, minFC) {
   df <- df[ c( groupping[[gr1]], groupping[[gr2]]  ) ]
   #Log2 fold change group2 - group1
   df$Log2FC <- apply(
     df, 1, function(x) {
       mean( x[ groupping[[gr2]] ] ) - mean( x[ groupping[[gr1]] ] )
     }
   )
   
   #T-test with equal variance
   df$T_Pval <- apply(
     df, 1, function(x) {
       res <- t.test(
         x[ groupping[[gr2]] ], x[ groupping[[gr1]] ],
         alternative = "two.sided", var.equal = TRUE
       )
       mean( x[ groupping[[gr2]] ] ) - mean( x[ groupping[[gr1]] ] )
       res$p.value
     }
   )
   #Benjamini-Hochberg correction for multiple testing
   df$adjPval <- p.adjust(df$T_Pval, method = "BH")
   df$Log10adjPval <- -1*log10(df$adjPval)
   #Add the categorical column for easier visualization
   df$Diff_Abund <- apply(
     df, 1, function(x) {
       if (x[["adjPval"]] <= maxAdjP & x[["Log2FC"]] >= minFC) {
         return( paste("Up in", gr2) )
       } else if (x[["adjPval"]] <= maxAdjP & x[["Log2FC"]] <= -1*minFC) {
         return( paste("Up in", gr1) )
       } else {
         return('Non-significant')
       }
     }
   )
   df
 }
maxAdjP <- 0.05
minLog2FC <- round(log2(1.3), 3)
gr1 <- "C"
gr2 <- "H"
dfTtest <- calc_ttest(countdata, groups, gr1, gr2, maxAdjP, minLog2FC)
 

#limma test of significant protein 
#P value calculated from https://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html
eb.fit <- function(dat, design){
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  logFC <- fit.eb$coefficients[, 2]
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, n)
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coefficients[, 2]/fit.eb$sigma/fit.eb$stdev.unscaled[, 2]
  t.mod <- fit.eb$t[, 2]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, 2]
  q.ord <- qvalue(p.ord)$q
  q.mod <- qvalue(p.mod)$q
  results.eb <- data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
  results.eb <- results.eb[order(results.eb$p.mod), ]
  return(results.eb)
}

tr <- c("H1", "H2", "H3")
ct <- c("C1", "C2", "C3")
design <- model.matrix(~factor(c(2,2,2,1,1,1)))
colnames(design) <- c("Intercept", "Diff")
res.eb <- eb.fit(countdata[, c(tr,ct)], design)


## Annoate the UniProt ID 
mouseUp <- UniProt.ws(10090)
annot <- select(
  x = mouseUp,
  columns <- c( "gene_primary"),
  keys = rownames(countdata),
  keytype = "UniProtKB"
)

annot <- annot[,-1] %>%
         rename(Accession=Entry, GeneID=Gene.Names..primary.) 

resannot <- as.data.frame(res.eb) %>% 
  rownames_to_column("Accession") %>% 
  left_join(annot, "Accession")

##Volcanoplot
png(file="Volcanoplot_peptide.png", width=800, height=800,res = 125)
ggplot(data = resannot, aes(x = logFC, y = -log10(p.mod))) +
  geom_point(data = subset(resannot), aes(size = abs(logFC)), color = "black", alpha = 0.1) +
  geom_point(data = subset(resannot, resannot$p.mod < 0.01 & resannot$logFC > 0.5), aes(size = abs(logFC)), color = "#ff0000", alpha = 0.3) +
  geom_point(data = subset(resannot, resannot$p.mod < 0.01 & resannot$logFC < -0.5), aes(size = abs(logFC)), color = "#56B1F7", alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05),lty = 4, lwd = 0.6, alpha = 0.8) +
  geom_vline(xintercept = c(0.5, -0.5), lty = 4, lwd = 0.6, alpha = 0.8) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  labs(x = "log2 (fold change)", y = "-log10 (moderated p-value)") +
  theme(legend.position = 'none') +
  geom_text_repel(data = rbind(subset(resannot, resannot$p.mod<0.01 & resannot$logFC > 0.5),
                               subset(resannot, resannot$p.mod<0.01 & resannot$logFC < -0.5)), 
                         aes(label = GeneID), 
                         col = "black", 
                         alpha = 0.8)	
dev.off()
                
##downstream functional analysis 
resannotsort <- resannot[order(resannot$p.mod),]
write.table(resannotsort, file="Peptide_DE.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

resannotsort_DN<-resannotsort %>% dplyr::filter(logFC <0.25 & p.mod < 0.01)
resannotsort_UP<-resannotsort %>% dplyr::filter(logFC >0.25 & p.mod < 0.01)

gene1_dn <- as.character(resannotsort_DN$GeneID)
# GO enrichment analysis
ego_BP_dn <- enrichGO(gene = gene1_dn,
                 OrgDb = org.Mm.eg.db, 
                 ont = "bp", 
                 keyType = 'SYMBOL',
                 pAdjustMethod = "fdr", 
                 pvalueCutoff = 0.05, 
                 minGSSize = 50, 
                 maxGSSize = 1750,
                 qvalueCutoff = 0.05,
                 readable = TRUE)

ego_BP_dn_simplify_table <- simplify(ego_BP_dn, cutoff=0.7, by="p.adjust", select_fun=min) %>%
                            as.data.frame()


ego_CC_dn <- enrichGO(gene = gene1_dn,
                      OrgDb = org.Mm.eg.db, 
                      ont = "cc", 
                      keyType = 'SYMBOL',
                      pAdjustMethod = "fdr", 
                      pvalueCutoff = 0.05, 
                      minGSSize = 50, 
                      maxGSSize = 1750,
                      qvalueCutoff = 0.05,
                      readable = TRUE)

ego_CC_dn_tree <- pairwise_termsim(ego_CC_dn)

png(file="CC_GO_treeplot_peptide.png", width=1800, height=1200, res = 150)
treeplot(ego_CC_dn_tree, showCategory = 25, color = "p.adjust", label_format = 30) +
  scale_color_gradient(limits = c(0, 0.05),low = "#56B1F7", high = "#228599") +
  theme(legend.position = "bottom")
dev.off()


gene1_up <- as.character(resannotsort_UP$GeneID)

ego_BP_up <- enrichGO(gene = gene1_up,
                      OrgDb = org.Mm.eg.db, 
                      ont = "bp", 
                      keyType = 'SYMBOL',
                      pAdjustMethod = "fdr", 
                      pvalueCutoff = 0.05, 
                      minGSSize = 50, 
                      maxGSSize = 1750,
                      qvalueCutoff = 0.05,
                      readable = TRUE)

ego_BP_up_simplify_table <- simplify(ego_BP_up, cutoff=0.7, by="p.adjust", select_fun=min) %>%
  as.data.frame()


#make a dataframe with the selected rows 
GoID_dn = c('GO:0006839',
            'GO:0046034',
            'GO:0046907',
            'GO:0070727',
            'GO:0007005')

selected_row <- data.frame()
for(i in GoID_dn){
  tmp  <- ego_BP_dn_simplify_table[ego_BP_dn_simplify_table$ID == i,]
  selected_row <- rbind(selected_row, tmp)
  }

png(file="dotplot_peptide_downregulated.png", width=1300, height=700, res = 120)
ggplot(selected_row,
  aes(x = GeneRatio, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = pvalue)) +
  scale_color_gradient(low = "#56B1F7", high = "#132B43") +
  ylab(NULL) +
  theme_bw() +
  ggtitle("Downgulated Pathways") +
  theme(axis.text.x = element_text(size=15,  color = "black"),
        axis.text.y = element_text(size=15, face="bold", color = "black"),
        plot.title = element_text(size=15, face = "bold")) 
dev.off()



png(file="dotplot_peptide.png", width=2000, height=700, res = 120)
cowplot::plot_grid(Downregulated_pathway, Upregulated_pathway)
dev.off()

# GSEA 
resannotsort_GSEA <-resannotsort %>% dplyr::filter(p.mod < 0.01)
## feature 1: numeric vector
gene_all = resannotsort_GSEA$logFC
## feature 2: named vector
names(gene_all) = as.character(resannotsort_GSEA$GeneID)
## feature 3: decreasing orde
gene_all = sort(gene_all, decreasing = TRUE)

#GESA for DE genes 
hallmark <- read.gmt("mh.all.v2022.1.Mm.symbols.gmt.txt")
c2 <- read.gmt("m2.all.v2022.1.Mm.symbols.gmt")
c5 <- read.gmt("m5.all.v2022.1.Mm.symbols.gmt")

egmt_c5 <- GSEA(gene_all, TERM2GENE = c5, pvalueCutoff = 1)  

geneSetID_c5 = c('GOCC_MITOCHONDRIAL_MATRIX',
                       'GOCC_MITOCHONDRIAL_ENVELOPE',
                       'GOCC_MITOCHONDRION')

png(file="C5_mito_GSEA_peptide.png", width=1000, height=700, res=140)
gseaNb(object = egmt_c5, 
       geneSetID = geneSetID_c5, 
       subPlot = 2,
       rmHt = T,
       termWidth = 35,
       addPval = T)
dev.off()


