library(GEOquery)
library(limma)
library(tidyverse)
library(ggfortify)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
library(fgsea)
library(labeling)
library(Hmisc)
library(ggrepel)

gset <- getGEO("GSE54646", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL4685", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("40545005520666533255040560035542044300050053206030",
               "20600545360653346500555666365665665665555661111111",
               "1111XXXXXXXXXXXX")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("PV_JAK2_Positive","Health","PV_JAK2_Negative","MF_JAK2_Negative","MF_JAK2_Positive","ET_JAK2_Positive","ET_JAK2_Negative"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)


fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cont.matrix <- makeContrasts(PVvsCon = PV_JAK2_Positive-Health, 
                             ETvsCon = ET_JAK2_Positive-Health,
                             MFvsCon = MF_JAK2_Positive-Health,
                             PVvsCon2 = PV_JAK2_Negative-Health,
                             ETvsCon2 = ET_JAK2_Negative-Health,
                             MFvsCon2 = MF_JAK2_Negative-Health,
                             levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT_PVvsCon <- topTable(fit2, coef=c("PVvsCon"), adjust="fdr", sort.by="logFC", p.value = 0.1, n = Inf)
Table_PVvsCon <- subset(tT_PVvsCon, select=c("Gene.symbol","adj.P.Val","logFC"))
resannot_PV <- Table_PVvsCon[!duplicated(Table_PVvsCon$Gene.symbol), ]
resannot_PV<-resannot_PV %>% dplyr::filter(!is.na(adj.P.Val) & !is.na(Gene.symbol))
# write.table(tT, file=stdout(), row.names=F, sep="\t")
trend_PV <-sapply(resannot_PV$logFC, function(x){if(x>0) 'up' else 'down'})
resannot_PV <- resannot_PV[order(resannot_PV$logFC, decreasing = T),]
pre_ranked_sig_genes_PV <- data.frame(resannot_PV, 'trend' = trend_PV, 'rank' = 1:nrow(resannot_PV), stringsAsFactors = F)
to_be_point_out_PV <- rbind(pre_ranked_sig_genes_PV[pre_ranked_sig_genes_PV$Gene.symbol == "CD24", ], 
                            pre_ranked_sig_genes_PV[pre_ranked_sig_genes_PV$Gene.symbol == "CD47", ])

tiff(file="fallplot_PV.tiff", width=700, height=700, res = 150)
  ggplot(pre_ranked_sig_genes_PV, aes(x=rank, y=logFC, color=logFC)) +
    geom_point(size=1)+
    geom_hline(yintercept = c(2,-2), linetype=2, size=0.25)+
    geom_hline(yintercept = c(0), linetype=1, size=0.5)+
    geom_vline(xintercept = 1636.5, linetype=2, size=0.25)+
    scale_color_gradient2(low="navy", high="firebrick3", mid="white", midpoint = 0)+
    geom_point(inherit.aes = F, data=to_be_point_out_PV, aes(x=rank, y=logFC), size = 3, color = 'black')+
    geom_point(inherit.aes = F, data=to_be_point_out_PV, aes(x=rank, y=logFC), size = 2, color = 'yellow')+
    ggrepel::geom_text_repel(inherit.aes = F, data = to_be_point_out_PV, aes(x=rank, y=logFC, label=Gene.symbol), size =5)+
    xlab('rank of differentially expressed genes') +
    theme_bw()+
    theme(panel.grid = element_line(color = 'white'), legend.title.align = 0.5)
  dev.off()
  
tT_ETvsCon <- topTable(fit2, coef=c("ETvsCon"), adjust="fdr", sort.by="logFC", p.value = 0.1, n = Inf)
Table_ETvsCon <- subset(tT_ETvsCon, select=c("Gene.symbol","adj.P.Val","logFC"))
resannot_ET <- Table_ETvsCon[!duplicated(Table_ETvsCon$Gene.symbol), ]
resannot_ET <-resannot_ET %>% dplyr::filter(!is.na(adj.P.Val) & !is.na(Gene.symbol))
trend_ET <-sapply(resannot_ET$logFC, function(x){if(x>0) 'up' else 'down'})
resannot_ET <- resannot_ET[order(resannot_ET$logFC, decreasing = T),]
pre_ranked_sig_genes_ET <- data.frame(resannot_ET, 'trend' = trend_ET, 'rank' = 1:nrow(resannot_ET), stringsAsFactors = F)
to_be_point_out_ET <- rbind(pre_ranked_sig_genes_ET[pre_ranked_sig_genes_ET$Gene.symbol == "CD24", ], pre_ranked_sig_genes_ET[pre_ranked_sig_genes_ET$Gene.symbol == "CD47", ])

tiff(file="fallplot_ET.tiff", width=700, height=700, res = 150)
ggplot(pre_ranked_sig_genes_ET, aes(x=rank, y=logFC, color=logFC)) +
  geom_point(size=1)+
  geom_hline(yintercept = c(2,-2), linetype=2, size=0.25)+
  geom_hline(yintercept = c(0), linetype=1, size=0.5)+
  geom_vline(xintercept = 1636.5, linetype=2, size=0.25)+
  scale_color_gradient2(low="navy", high="firebrick3", mid="white", midpoint = 0)+
  geom_point(inherit.aes = F, data=to_be_point_out_ET, aes(x=rank, y=logFC), size = 3, color = 'black')+
  geom_point(inherit.aes = F, data=to_be_point_out_ET, aes(x=rank, y=logFC), size = 2, color = 'yellow')+
  ggrepel::geom_text_repel(inherit.aes = F, data = to_be_point_out_ET, aes(x=rank, y=logFC, label=Gene.symbol), size =5)+
  xlab('rank of differentially expressed genes') +
  theme_bw()+
  theme(panel.grid = element_line(color = 'white'), legend.title.align = 0.5)
dev.off()

tT_MFvsCon <- topTable(fit2, coef=c("MFvsCon"), adjust="fdr", sort.by="logFC", p.value = 0.1, n = Inf)
Table_MFvsCon <- subset(tT_MFvsCon, select=c("Gene.symbol","adj.P.Val","logFC"))
resannot_MF <- Table_MFvsCon[!duplicated(Table_MFvsCon$Gene.symbol), ]
resannot_MF <-resannot_MF %>% dplyr::filter(!is.na(adj.P.Val) & !is.na(Gene.symbol))
trend_MF <-sapply(resannot_MF$logFC, function(x){if(x>0) 'up' else 'down'})
resannot_MF <- resannot_MF[order(resannot_MF$logFC, decreasing = T),]
pre_ranked_sig_genes_MF <- data.frame(resannot_MF, 'trend' = trend_MF, 'rank' = 1:nrow(resannot_MF), stringsAsFactors = F)
to_be_point_out_MF <- rbind(pre_ranked_sig_genes_MF[pre_ranked_sig_genes_MF$Gene.symbol == "CD24", ], pre_ranked_sig_genes_MF[pre_ranked_sig_genes_MF$Gene.symbol == "CD47", ])

tiff(file="fallplot_MF.tiff", width=700, height=700, res = 150)
ggplot(pre_ranked_sig_genes_MF, aes(x=rank, y=logFC, color=logFC)) +
  geom_point(size=1)+
  geom_hline(yintercept = c(2,-2), linetype=2, size=0.25)+
  geom_hline(yintercept = c(0), linetype=1, size=0.5)+
  geom_vline(xintercept = 1636.5, linetype=2, size=0.25)+
  scale_color_gradient2(low="navy", high="firebrick3", mid="white", midpoint = 0)+
  geom_point(inherit.aes = F, data=to_be_point_out_MF, aes(x=rank, y=logFC), size = 3, color = 'black')+
  geom_point(inherit.aes = F, data=to_be_point_out_MF, aes(x=rank, y=logFC), size = 2, color = 'yellow')+
  ggrepel::geom_text_repel(inherit.aes = F, data = to_be_point_out_MF, aes(x=rank, y=logFC, label=Gene.symbol), size =5)+
  xlab('rank of differentially expressed genes') +
  theme_bw()+
  theme(panel.grid = element_line(color = 'white'), legend.title.align = 0.5)
dev.off()


fit2 <- eBayes(fit2, 0.01)
tT_PVvsCon2 <- topTable(fit2, coef=c("PVvsCon2"), adjust="fdr", sort.by="logFC", p.value = 0.1, n = Inf)
Table_PVvsCon2 <- subset(tT_PVvsCon2, select=c("Gene.symbol","adj.P.Val","logFC"))
tT_ETvsCon2 <- topTable(fit2, coef=c("ETvsCon2"), adjust="fdr", sort.by="logFC", p.value = 0.1, n = Inf)
Table_ETvsCon2 <- subset(tT_ETvsCon2, select=c("Gene.symbol","adj.P.Val","logFC"))
tT_MFvsCon2 <- topTable(fit2, coef=c("MFvsCon2"), adjust="fdr", sort.by="logFC", p.value = 0.1, n = Inf)
Table_MFvsCon2 <- subset(tT_MFvsCon2, select=c("Gene.symbol","adj.P.Val","logFC"))