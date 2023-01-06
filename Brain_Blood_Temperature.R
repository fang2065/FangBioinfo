library(tidyverse)

## import data
sampleinfo <- read_tsv("Sampleinfo.txt")
seqdata <- read_tsv("Counts.txt", comment="#")
seqdata <- seqdata[!duplicated(seqdata$GeneID), ]

## reformat data
countdata <- seqdata %>%
  as.data.frame() %>% 
  column_to_rownames("GeneID") %>% # turn the geneid column into rownames
  dplyr::select(sampleinfo$Sample) %>% # keep sample columns using sampleinfo$Sample
  as.matrix()

## remove non-expressed genes
keep <- rowSums(countdata) > 24
countdata <- countdata[keep,]
countdata <- countdata[apply(countdata, 1, function(x) sum( is.na(x) ))==0,]

library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library(ggfortify)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
library(labeling)
library(Hmisc)

 
## bar plots library size
librarySizes <- colSums(countdata)
col <- brewer.pal(12, "Paired")
png(file="BarPlot.png", width=2000, height=1500, res = 200)
par(mar=c(6, 6, 4, 4), las=2)
par(mfrow=c(1,2))
barplot(librarySizes[1:12], 
        names=names(librarySizes[1:12]), 
        col=col)
abline(h=5.0e6, lty=2)
title(main="A. Blood")
barplot(librarySizes[13:24], 
        names=names(librarySizes[13:24]), 
        col=col)
abline(h=5.0e6, lty=2)
title(main="B. Brain")
dev.off()

lcpm_countdata <- cpm(countdata, log=TRUE)

nsamples <- ncol(countdata)
col <- brewer.pal(12, "Paired")
png(file="reads_density_plot.png", width=1600, height=1200, res = 200)
par(mfrow=c(1,2))
lcpm_countdata1 <- lcpm_countdata[,1:12]
plot(density(lcpm_countdata1[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Blood", xlab="Log-count")
for (i in 2:nsamples){
  den <- density(lcpm_countdata1[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(countdata[,1:12]), text.col=col,cex=0.5, bty="n")
lcpm_countdata2 <- lcpm_countdata[,13:24]
plot(density(lcpm_countdata2[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Brain", xlab="Log-count")
for (i in 2:nsamples){
  den <- density(lcpm_countdata2[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(countdata[,13:24]), text.col=col,cex=0.5, bty="n")
dev.off()

logcounts <- log2(countdata + 1)
# Check distributions of samples using boxplots
png(file="BoxCounts.png", width=2400, height=1200)
par(mar=c(6, 6, 4, 4), las=2)
par(mfrow=c(1,2))
boxplot(logcounts[,1:12], 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=col)
# Let's add a blue horizontal line that corresponds to the median logcounts
abline(h=median(as.matrix(logcounts)), col="blue")
title(main="A. Blood")
boxplot(logcounts[,13:24], 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=col)
# Let's add a blue horizontal line that corresponds to the median logcounts
abline(h=median(as.matrix(logcounts)), col="blue")
title(main="B. Brain")
dev.off()


## make a heatmap of correlation matrix now
cormat <- round(cor(logcounts, method = c("pearson")),2)
res <- rcorr(as.matrix(cormat)) #calculate the p value for each correlation coefficient
res$p #check the value  
ha_column = HeatmapAnnotation(df = data.frame(status = c(rep("Blood", 12), rep("Brain", 12))),
                              col = list(status = c("Blood" =  "orange", "Brain" = "purple")))
png(file="Correlation.png", width=1300, height=1500, res = 200)
par(mar=c(6, 6, 4, 4), las=2)
Heatmap(cormat, name = "PCC", column_title = "Pearson Correlation Matrix", top_annotation = ha_column)
dev.off()

#filer the data using filterByExpr in limma 
group <- as.factor(sampleinfo$Group)
keep.exprs <- filterByExpr(countdata, group = group)
countdata_filtered <- countdata[keep.exprs,]
logcountdata_filtered <- log2(countdata_filtered + 1)
lcpm_filtered <- cpm(countdata_filtered, log=TRUE)

#run qc again after filtering 
nsamples <- ncol(countdata)
col <- brewer.pal(12, "Paired")
png(file="reads_density_plot_filtered.png", width=1600, height=1200, res = 200)
par(mfrow=c(1,2))
lcpm_filtered1 <- lcpm_filtered[,1:12]
plot(density(lcpm_filtered1[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Blood", xlab="Lcpm")
for (i in 2:nsamples){
  den <- density(lcpm_filtered1[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(countdata[,1:12]), text.col=col,cex=0.5, bty="n")
lcpm_filtered2 <- lcpm_filtered[,13:24]
plot(density(lcpm_filtered2[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Brain", xlab="Lcpm")
for (i in 2:nsamples){
  den <- density(lcpm_filtered2[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(countdata[,13:24]), text.col=col,cex=0.5, bty="n")
dev.off()

png(file="BoxCounts.png", width=2400, height=1200)
par(mar=c(6, 6, 4, 4), las=2)
par(mfrow=c(1,2))
boxplot(logcountdata_filtered[,1:12], 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=col)
abline(h=median(as.matrix(logcounts)), col="blue")
title(main="A. Blood")
boxplot(logcountdata_filtered[,13:24], 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=col)
abline(h=median(as.matrix(logcounts)), col="blue")
title(main="B. Brain")
dev.off()



cormat_filtered <- round(cor(logcountdata_filtered, method = c("pearson")),2)
ha_column = HeatmapAnnotation(df = data.frame(status = c(rep("Blood", 6), rep("Brain", 6))),
                              col = list(status = c("Blood" =  "orange", "Brain" = "purple")))
png(file="Correlation_filtered.png", width=1300, height=1500, res = 200)
par(mar=c(6, 6, 4, 4), las=2)
Heatmap(cormat_filtered, name = "PCC", column_title = "Pearson Correlation Matrix", top_annotation = ha_column)
dev.off()

pcDat <- prcomp(t(logcountdata_filtered))
# plot PCA
png(file="PCA12_tissue.png", width=1200, height=1200)
autoplot(pcDat, data = sampleinfo, fill="Tissue", shape=21, size=5) + scale_shape_manual(values=c(21)) + guides(fill = guide_legend(override.aes=list(shape=21)))
dev.off()
png(file="PCA23_tissue.png", width=1200, height=1200)
autoplot(pcDat, data = sampleinfo, fill="Tissue", shape=21, size=5, x=2, y=3) + scale_shape_manual(values=c(21)) + guides(fill = guide_legend(override.aes=list(shape=21)))
dev.off()
png(file="PCA13_tissue.png", width=1200, height=1200)
autoplot(pcDat, data = sampleinfo, fill="Tissue", shape=21, size=5, x=1, y=3) + scale_shape_manual(values=c(21)) + guides(fill = guide_legend(override.aes=list(shape=21)))
dev.off()

png(file="PCA12_Group.png", width=1200, height=1200)
autoplot(pcDat, data = sampleinfo, fill="Group", shape=21, size=5) + scale_shape_manual(values=c(21)) + guides(fill = guide_legend(override.aes=list(shape=21)))
dev.off()
png(file="PCA23_Group.png", width=1200, height=1200)
autoplot(pcDat, data = sampleinfo, fill="Group", shape=21, size=5, x=2, y=3) + scale_shape_manual(values=c(21)) + guides(fill = guide_legend(override.aes=list(shape=21)))
dev.off()
png(file="PCA13_Group.png", width=1200, height=1200)
autoplot(pcDat, data = sampleinfo, fill="Group", shape=21, size=5, x=1, y=3) + scale_shape_manual(values=c(21)) + guides(fill = guide_legend(override.aes=list(shape=21)))
dev.off()


# Let's analyse blood samples first...
countdata_blood <- countdata_filtered[,1:12]
sampleinfo_blood <- sampleinfo[1:12,]
group <- as.factor(sampleinfo_blood$Group)
design <- model.matrix(~ 0 + group)
colnames(design) <- gsub("group", "", colnames(design))
contr.matrix <- makeContrasts(
  ConTvsHighT = C - H, 
  ConTvsHighT_E = C - HE, 
  HighTvsHighT_E = H - HE, 
  levels = colnames(design))

png(file="blood_voom_adjustment.png", width=2500, height=1500, res = 200)
par(mfrow=c(1,2))
countdata_blood_voom <- voom(countdata_blood, design, plot=TRUE)
countdata_blood_voom_fit <- lmFit(countdata_blood_voom, design)
countdata_blood_voom_fit <- contrasts.fit(countdata_blood_voom_fit, contrasts=contr.matrix)
countdata_blood_voom_efit <- eBayes(countdata_blood_voom_fit)
plotSA(countdata_blood_voom_efit, main="Final model: Mean-variance trend")
dev.off()

summary(decideTests(countdata_blood_voom_efit))
countdata_blood_voom_tfit <- treat(countdata_blood_voom_fit, lfc=1)
dt_blood <- decideTests(countdata_blood_voom_tfit)
circleCol <- brewer.pal(3, "Pastel2")
png(file="venn_blood_limma.png", width=1500, height=1500, res = 200)
vennDiagram(dt[,1:3], circle.col=circleCol)
dev.off()


# ...and then brain samples

countdata_brain <- countdata_filtered[,13:24]
sampleinfo_brain <- sampleinfo[13:24,]
group <- as.factor(sampleinfo_brain$Group)
design <- model.matrix(~ 0 + group)
colnames(design) <- gsub("group", "", colnames(design))
contr.matrix <- makeContrasts(
  ConTvsHighT = C - H, 
  ConTvsHighT_E = C - HE, 
  HighTvsHighT_E = H - HE, 
  levels = colnames(design))

png(file="brain_voom_adjustment.png", width=2500, height=1500, res = 200)
par(mfrow=c(1,2))
countdata_brain_voom <- voom(countdata_brain, design, plot=TRUE)
countdata_brain_voom_fit <- lmFit(countdata_brain_voom, design)
countdata_brain_voom_fit <- contrasts.fit(countdata_brain_voom_fit, contrasts=contr.matrix)
countdata_brain_voom_efit <- eBayes(countdata_brain_voom_fit)
plotSA(countdata_brain_voom_efit, main="Final model: Mean-variance trend")
dev.off()

summary(decideTests(countdata_brain_voom_efit))
countdata_brain_voom_tfit <- treat(countdata_brain_voom_fit, lfc=1)
dt_brain <- decideTests(countdata_brain_voom_tfit)
circleCol <- brewer.pal(3, "Pastel2")
png(file="venn_brain_limma.png", width=1500, height=1500, res = 200)
vennDiagram(dt_brain[,1:3], circle.col=circleCol)
dev.off()

write.fit(countdata_brain_voom_tfit, dt_brain, file="brain_limma_results.txt")
ConTvsHighT <- topTreat(countdata_brain_voom_tfit, coef=1, n=Inf)
ConTvsHighT_E <- topTreat(countdata_brain_voom_tfit, coef=2, n=Inf)
HighTvsHighT_E <- topTreat(countdata_brain_voom_tfit, coef=3, n=Inf)

#Annoate 
ensembl=useMart("ENSEMBL_MART_ENSEMBL", host="useast.ensembl.org")
ensembl = useDataset("mmusculus_gene_ensembl", mart=ensembl)
attributeNames <- c('ensembl_gene_id', 'entrezgene_id', 'mgi_symbol', 'description', 'gene_biotype', 'chromosome_name', 'start_position', 'end_position', 'strand')
ourFilterType <- "ensembl_gene_id"
filterValues <- rownames(ConTvsHighT_E)
annot <- getBM(attributes=attributeNames, filters = ourFilterType, values = filterValues, mart = ensembl)
names(annot) <- c("GeneID", "Entrez", "Symbol", "Description", "Biotype", "Chr", "Start", "End", "Strand")
resannot <- as.data.frame(res) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(annot, "GeneID")

#DESeq2_blood
library(DESeq2)
design <- as.formula(~ 0 + Group)
ddsObj.raw_blood <- DESeqDataSetFromMatrix(countData = countdata_blood, colData = sampleinfo_blood, design = design)
ddsObj_blood <- DESeq(ddsObj.raw_blood)
resultsNames(ddsObj_blood)
res_CvsH_blood <- results(ddsObj_blood, contrast=c("Group", "C", "H"), alpha=0.05)
res_CvsHE_blood <- results(ddsObj_blood, contrast=c("Group", "C", "HE"), alpha=0.05)
res_HvsHE_blood <- results(ddsObj_blood, contrast=c("Group", "H", "HE"), alpha=0.05)

#Annoate the genelist 
ensembl=useMart("ENSEMBL_MART_ENSEMBL", host="https://useast.ensembl.org")
ensembl = useDataset("mmusculus_gene_ensembl", mart=ensembl)
attributeNames <- c('ensembl_gene_id', 'entrezgene_id', 'mgi_symbol', 'description', 'gene_biotype', 'chromosome_name', 'start_position', 'end_position', 'strand')
ourFilterType <- "ensembl_gene_id"
filterValues <- rownames(res_CvsH_blood)
annot <- getBM(attributes=attributeNames, filters = ourFilterType, values = filterValues, mart = ensembl)
names(annot) <- c("GeneID", "Entrez", "Symbol", "Description", "Biotype", "Chr", "Start", "End", "Strand")
resannot_CvsH_blood <- as.data.frame(res_CvsH_blood) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(annot, "GeneID")

ensembl=useMart("ENSEMBL_MART_ENSEMBL", host="https://useast.ensembl.org")
ensembl = useDataset("mmusculus_gene_ensembl", mart=ensembl)
attributeNames <- c('ensembl_gene_id', 'entrezgene_id', 'mgi_symbol', 'description', 'gene_biotype', 'chromosome_name', 'start_position', 'end_position', 'strand')
ourFilterType <- "ensembl_gene_id"
filterValues <- rownames(res_CvsHE_blood)
annot <- getBM(attributes=attributeNames, filters = ourFilterType, values = filterValues, mart = ensembl)
names(annot) <- c("GeneID", "Entrez", "Symbol", "Description", "Biotype", "Chr", "Start", "End", "Strand")
resannot_CvsHE_blood <- as.data.frame(res_CvsHE_blood) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(annot, "GeneID")

ensembl=useMart("ENSEMBL_MART_ENSEMBL", host="https://useast.ensembl.org")
ensembl = useDataset("mmusculus_gene_ensembl", mart=ensembl)
attributeNames <- c('ensembl_gene_id', 'entrezgene_id', 'mgi_symbol', 'description', 'gene_biotype', 'chromosome_name', 'start_position', 'end_position', 'strand')
ourFilterType <- "ensembl_gene_id"
filterValues <- rownames(res_HvsHE_blood)
annot <- getBM(attributes=attributeNames, filters = ourFilterType, values = filterValues, mart = ensembl)
names(annot) <- c("GeneID", "Entrez", "Symbol", "Description", "Biotype", "Chr", "Start", "End", "Strand")
resannot_HvsHE_blood <- as.data.frame(res_HvsHE_blood) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(annot, "GeneID")

## export results
resannotsort_CvsH_blood <- resannot_CvsH_blood[order(resannot_CvsH_blood$padj),] %>% 
  dplyr::filter(Biotype == "protein_coding") %>%
  dplyr::filter(!is.na(padj) & !is.na(Entrez))

resannotsort_CvsHE_blood <- resannot_CvsHE_blood[order(resannot_CvsHE_blood$padj),]%>%
  dplyr::filter(Biotype == "protein_coding") %>%
  dplyr::filter(!is.na(padj) & !is.na(Entrez))
  
resannotsort_HvsHE_blood <- resannot_HvsHE_blood[order(resannot_HvsHE_blood$padj),]%>%
  dplyr::filter(Biotype == "protein_coding") %>%
  dplyr::filter(!is.na(padj) & !is.na(Entrez))
  
  
write.table(resannotsort_CvsH_blood, file="Gene_DE_CvsH_blood.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(resannotsort_CvsHE_blood, file="Gene_DE_CvsHE_blood.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(resannotsort_HvsHE_blood, file="Gene_DE_HvsHE_blood.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


ddsObj.raw_brain <- DESeqDataSetFromMatrix(countData = countdata_brain, colData = sampleinfo_brain, design = design)
ddsObj_brain <- DESeq(ddsObj.raw_brain)
resultsNames(ddsObj_brain)
res_CvsH_brain <- results(ddsObj_brain, contrast=c("Group", "C", "H"), alpha=0.05)
res_CvsHE_brain <- results(ddsObj_brain, contrast=c("Group", "C", "HE"), alpha=0.05)
res_HvsHE_brain <- results(ddsObj_brain, contrast=c("Group", "H", "HE"), alpha=0.05)

ensembl=useMart("ENSEMBL_MART_ENSEMBL", host="https://useast.ensembl.org")
ensembl = useDataset("mmusculus_gene_ensembl", mart=ensembl)
attributeNames <- c('ensembl_gene_id', 'entrezgene_id', 'mgi_symbol', 'description', 'gene_biotype', 'chromosome_name', 'start_position', 'end_position', 'strand')
ourFilterType <- "ensembl_gene_id"
filterValues <- rownames(res_CvsH_brain)
annot <- getBM(attributes=attributeNames, filters = ourFilterType, values = filterValues, mart = ensembl)
names(annot) <- c("GeneID", "Entrez", "Symbol", "Description", "Biotype", "Chr", "Start", "End", "Strand")
resannot_CvsH_brain <- as.data.frame(res_CvsH_brain) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(annot, "GeneID")

ensembl=useMart("ENSEMBL_MART_ENSEMBL", host="https://useast.ensembl.org")
ensembl = useDataset("mmusculus_gene_ensembl", mart=ensembl)
attributeNames <- c('ensembl_gene_id', 'entrezgene_id', 'mgi_symbol', 'description', 'gene_biotype', 'chromosome_name', 'start_position', 'end_position', 'strand')
ourFilterType <- "ensembl_gene_id"
filterValues <- rownames(res_CvsHE_brain)
annot <- getBM(attributes=attributeNames, filters = ourFilterType, values = filterValues, mart = ensembl)
names(annot) <- c("GeneID", "Entrez", "Symbol", "Description", "Biotype", "Chr", "Start", "End", "Strand")
resannot_CvsHE_brain <- as.data.frame(res_CvsHE_brain) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(annot, "GeneID")

ensembl=useMart("ENSEMBL_MART_ENSEMBL", host="https://useast.ensembl.org")
ensembl = useDataset("mmusculus_gene_ensembl", mart=ensembl)
attributeNames <- c('ensembl_gene_id', 'entrezgene_id', 'mgi_symbol', 'description', 'gene_biotype', 'chromosome_name', 'start_position', 'end_position', 'strand')
ourFilterType <- "ensembl_gene_id"
filterValues <- rownames(res_HvsHE_brain)
annot <- getBM(attributes=attributeNames, filters = ourFilterType, values = filterValues, mart = ensembl)
names(annot) <- c("GeneID", "Entrez", "Symbol", "Description", "Biotype", "Chr", "Start", "End", "Strand")
resannot_HvsHE_brain <- as.data.frame(res_HvsHE_brain) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(annot, "GeneID")

resannotsort_CvsH_brain <- resannot_CvsH_brain[order(resannot_CvsH_brain$padj),] %>% 
  dplyr::filter(Biotype == "protein_coding") %>%
  dplyr::filter(!is.na(padj) & !is.na(Entrez))

resannotsort_CvsHE_brain <- resannot_CvsHE_brain[order(resannot_CvsHE_brain$padj),]%>%
  dplyr::filter(Biotype == "protein_coding") %>%
  dplyr::filter(!is.na(padj) & !is.na(Entrez))

resannotsort_HvsHE_brain <- resannot_HvsHE_brain[order(resannot_HvsHE_brain$padj),]%>%
  dplyr::filter(Biotype == "protein_coding") %>%
  dplyr::filter(!is.na(padj) & !is.na(Entrez))


write.table(resannotsort_CvsH_brain, file="Gene_DE_CvsH_brain.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(resannotsort_CvsHE_brain, file="Gene_DE_CvsHE_brain.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(resannotsort_HvsHE_brain, file="Gene_DE_HvsHE_brain.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

countdata_blood_normalised <- counts(ddsObj_blood, normalized=TRUE)
write.table(countdata_blood_normalised, file="Count_blood_normalised.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

countdata_brain_normalised <- counts(ddsObj_brain, normalized=TRUE)
write.table(countdata_brain_normalised, file="Count_brain_normalised.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
