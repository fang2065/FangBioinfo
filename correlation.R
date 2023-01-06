library(DaMiRseq)
library(tidyverse)
library(DESeq2)
library(glmnet)
library(edgeR)
library(mixOmics)
library(WGCNA)
## import data
seqdata <- read_tsv("counts_PMF.txt", comment="#")
seqdata <- seqdata[!duplicated(seqdata$Geneid), ]

countdata <- seqdata %>%
  as.data.frame() %>% 
  column_to_rownames("Geneid") %>% # turn the geneid column into rownames
  as.matrix()

sampleinfo_update <- read_tsv("EGA_Samples_PMF.txt")
covariatedata<- sampleinfo_update %>%
  as.data.frame() %>%  
  column_to_rownames("Patient") 

covariatedata$class <- as.numeric(covariatedata$class)
covariatedata$mutation <- as.factor(covariatedata$mutation)

keep.exprs <- filterByExpr(countdata)
countdata_filtered <- countdata[keep.exprs,]

design <- as.formula(~ mutation)
covariatedata$mutation <- factor(covariatedata$mutation, levels = c(rev(unique(covariatedata$mutation))))
ddsObj.raw <- DESeqDataSetFromMatrix(countData = countdata_filtered, colData = covariatedata, design = design)
ddsObj.raw <- estimateSizeFactors(ddsObj.raw)
counts_normalised <- counts(ddsObj.raw, normalized=TRUE)


sampleTree = hclust(dist(t(counts_normalised)), method = "average")
png(file="hclust.png", width=1200, height=1200)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
abline(h=900000, col="blue")
dev.off()

clust = cutreeStatic(sampleTree, cutHeight = 900000, minSize = 10)
table(clust)
keepSamples = (clust==1)
counts_coll_filtered = counts_normalised[,keepSamples]

sampleTree = hclust(dist(t(counts_coll_filtered)), method = "average")
png(file="hclust_filtered.png", width=1200, height=1200)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
abline(h=900000, col="blue")
dev.off()

covariatedata_filter = covariatedata[colnames(counts_coll_filtered),]

x <- log10(counts_coll_filtered[rownames(counts_coll_filtered) == "ENSG00000272398",])
y <- covariatedata_filter$class

table <- t(rbind(x, y)) %>%
  as.data.frame()
colnames(table)[1]<-'log10(CD24 expression level)'
colnames(table)[2]<-'JAK2V617F mutation allele burden'

png(file="spearman correlation in PMF.png", width=1000, height=1000, res=150)
ggscatter(table, 
          x = "x",
          y = "y",
          add = "reg.line", 
          size = 3,
          shape = 21,
          color = "red",
          add.params = list(color = "red",fill = "lightgrey", size = 1),
          conf.int = TRUE, 
          cor.coef = TRUE, 
          cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n"),
          main="Primary Myelofibrosis",
          xlab="log10(CD24 expression level)",
          ylab="JAK2V617F mutation allele burden (%)")
dev.off()

MyResult.splsda <- splsda(X, Y, keepX = c(50,50)) 
plotIndiv(MyResult.splsda)
plotVar(MyResult.splsda)

plotIndiv(MyResult.splsda, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = TRUE, title = 'sPLS-DA',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')

