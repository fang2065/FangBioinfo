setwd("~/Documents/FW")

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


#filer the data using filterByExpr in limma 
group <- as.factor(sampleinfo$Group)
keep.exprs <- filterByExpr(countdata, group = group)
countdata_filtered <- countdata[keep.exprs,]
logcountdata_filtered <- log2(countdata_filtered + 1)
lcpm_filtered <- cpm(countdata_filtered, log=TRUE)

#DESeq2_vrain
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

library(DESeq2)
design <- as.formula(~ 0 + Group)
ddsObj.raw_brain <- DESeqDataSetFromMatrix(countData = countdata_brain, colData = sampleinfo_brain, design = design)
ddsObj_brain <- DESeq(ddsObj.raw_brain)
resultsNames(ddsObj_brain)
res_CvsH_brain <- results(ddsObj_brain, contrast=c("Group", "H", "C"), alpha=0.05)


#Annoate the genelist 
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


## export results
resannotsort_CvsH_brain <- resannot_CvsH_brain[order(resannot_CvsH_brain$padj),] %>% 
  dplyr::filter(Biotype == "protein_coding") %>%
  dplyr::filter(!is.na(padj) & !is.na(Entrez))
  

library(DOSE)
library(org.Mm.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(ggrepel)
library(ggridges)
library(GseaVis)



png(file="Volcanoplot.png", width=800, height=800,res = 125)
ggplot(data = resannotsort_CvsH_brain, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(data = subset(resannotsort_CvsH_brain), aes(size = abs(log2FoldChange)), color = "black", alpha = 0.1) +
  geom_point(data = subset(resannotsort_CvsH_brain, resannotsort_CvsH_brain$pvalue<0.05 & resannotsort_CvsH_brain$log2FoldChange > 2), aes(size = abs(log2FoldChange)), color = "#ff0000", alpha = 0.3) +
  geom_point(data = subset(resannotsort_CvsH_brain, resannotsort_CvsH_brain$pvalue<0.05 & resannotsort_CvsH_brain$log2FoldChange < -2), aes(size = abs(log2FoldChange)), color = "#56B1F7", alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05),lty = 4, lwd = 0.6, alpha = 0.8) +
  geom_vline(xintercept = c(1, -1), lty = 4, lwd = 0.6, alpha = 0.8) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  labs(x = "log2 (fold change)", y = "-log10 (pvalue)") +
  theme(legend.position = 'none') +
  geom_text_repel(data = rbind(subset(resannotsort_CvsH_brain, resannotsort_CvsH_brain$pvalue<0.05 & resannotsort_CvsH_brain$log2FoldChange > 2),
                               subset(resannotsort_CvsH_brain, resannotsort_CvsH_brain$pvalue<0.05 & resannotsort_CvsH_brain$log2FoldChange < -2)), 
                  aes(label = Symbol), 
                  col = "black", 
                  alpha = 0.8)	
dev.off()

resannotsort_DN<-resannotsort_CvsH_brain %>% dplyr::filter(log2FoldChange<0 & pvalue < 0.05)
resannotsort_UP<-resannotsort_CvsH_brain %>% dplyr::filter(log2FoldChange>0 & pvalue < 0.05)


gene1_dn <- as.character(resannotsort_DN$Entrez)

#for example 
ego_BP_dn <- enrichGO(gene = gene1_dn, 
                       OrgDb = org.Mm.eg.db, 
                       ont = "bp", 
                       pAdjustMethod = "fdr", 
                       pvalueCutoff = 0.05, 
                       minGSSize = 50, 
                       maxGSSize = 1750,
                       qvalueCutoff = 0.05,
                       readable = TRUE)


ego_BP_dn_simplify_table <- simplify(ego_BP_dn, cutoff=0.7, by="p.adjust", select_fun=min) %>%
                            as.data.frame()

Downregulated_pathway <- ggplot(ego_BP_dn_simplify_table[c(1,6,8,9,10),],
                         aes(x = GeneRatio, y = Description)) + 
                         geom_point(aes(size = GeneRatio, color = pvalue)) +
                         scale_color_gradient(low = "#56B1F7", high = "#132B43") +
                         ylab(NULL) +
                         theme_bw() +
                         ggtitle("Downgulated Pathways") +
                         theme(axis.text.x = element_text(size=15,  color = "black"),
                               axis.text.y = element_text(size=15, face="bold", color = "black"),
                               plot.title = element_text(size=15, face = "bold")) 
                         

gene1_up <- as.character(resannotsort_UP$Entrez)

ego_BP_up <- enrichGO(gene = gene1_up, 
                      OrgDb = org.Mm.eg.db, 
                      ont = "bp", 
                      pAdjustMethod = "fdr", 
                      pvalueCutoff = 0.05, 
                      minGSSize = 50, 
                      maxGSSize = 1750,
                      qvalueCutoff = 0.05,
                      readable = TRUE)

ego_BP_up_simplify_table <- simplify(ego_BP_up, cutoff=0.7, by="p.adjust", select_fun=min) %>%
                            as.data.frame() 

Upregulated_pathway <- ggplot(ego_BP_up_simplify_table[c(3,6,10,24,137),],
                              aes(x = GeneRatio, y = Description)) + 
                              geom_point(aes(size = GeneRatio , color = pvalue)) +
                              scale_size_area(max_size = 5) +
                              scale_color_gradient(low = "#ff3232", high = "#132B43") +
                              ylab(NULL) +
                              theme_bw() +
                              ggtitle("Upregulated Pathways") +
                              theme(axis.text.x = element_text(size=15,  color = "black"),
                                    axis.text.y = element_text(size=15, face="bold", color = "black"),
                                    plot.title = element_text(size=15, face = "bold")) 
                              


png(file="dotplot.png", width=2000, height=700, res = 120)
cowplot::plot_grid(Downregulated_pathway, Upregulated_pathway)
dev.off()



resannotsort_GSEA <-resannotsort_CvsH_brain %>% dplyr::filter(pvalue < 0.05)
## feature 1: numeric vector
gene_all = resannotsort_GSEA$log2FoldChange
## feature 2: named vector
names(gene_all) = as.character(resannotsort_GSEA$Entrez)
## feature 3: decreasing orde
gene_all = sort(gene_all, decreasing = TRUE)

#GESA for DE genes 
hallmark <- read.gmt("mh.all.v2022.1.Mm.entrez.gmt")
c2 <- read.gmt("m2.all.v2022.1.Mm.entrez.gmt")
c5 <- read.gmt("m5.all.v2022.1.Mm.entrez.gmt")
c8 <- read.gmt("m8.all.v2022.1.Mm.entrez.gmt")

egmt_c8_singlecell <- GSEA(gene_all, TERM2GENE = c8, pvalueCutoff = 1)  

png(file="ridgeplot.png", width=1200, height=1200)
ridgeplot(egmt_c8_singlecell, showCategory=7) + 
  scale_fill_gradient(low = "#ff3232", high = "#b22323") +
  xlab('NES') +
  theme(axis.text.x = element_text(size=15,  color = "black"),
        axis.text.y = element_text(size=15, face="bold", color = "black"),
        plot.title = element_text(size=15, face = "bold")) 
dev.off()

egmt_hallmark <- GSEA(gene_all, TERM2GENE = hallmark, pvalueCutoff = 1)  

geneSetID_hallmark = c('HALLMARK_TNFA_SIGNALING_VIA_NFKB',
              'HALLMARK_INFLAMMATORY_RESPONSE',
              'HALLMARK_IL2_STAT5_SIGNALING')

png(file="hallmark_GSEA.png", width=1000, height=700, res=140)
gseaNb(object = egmt_hallmark, 
       geneSetID = geneSetID_hallmark, 
       subPlot = 2,
       rmHt = T,
       termWidth = 35,
       addPval = T)
dev.off()

egmt_c2 <- GSEA(gene_all, TERM2GENE = c2, pvalueCutoff = 1)  
egmt_c5 <- GSEA(gene_all, TERM2GENE = c5, pvalueCutoff = 1)  

geneSetID_c2 = c('REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS')

png(file="hallmark_c2.png", width=700, height=700, res=140)
gseaNb(object = egmt_c2, 
       geneSetID = geneSetID_c2, 
       subPlot = 2,
       rmHt = T,
       termWidth = 35,
       addPval = T)
dev.off()

#raw a heat map using a specific gene set

ap1_raw <- read_tsv("mito.txt")
ap1 <- t(ap1_raw$converted_alias) #converted_alias here is EntrezGeneIds
symbols <- as.character(c(ap1))
Nsig1 <- length(symbols)

plotDat <- vst(ddsObj_brain)[symbols,] %>% assay()

#transfer EntrezGeneIds into Gene symbol 
rownames(plotDat) <- c(ap1_raw$initial_alias)

z.mat <- t(scale(t(plotDat[,1:6]), center=TRUE, scale=TRUE))

myPalette <- c("blue", "white", "red")
myRamp = colorRamp2(c(-2, 0, 2), myPalette)

hcDat <- hclust(dist(z.mat))
cutGroups <- cutree(hcDat, h=4)
ha1 = ha_column = HeatmapAnnotation(df = data.frame(status = c(rep("Control", 3), rep("High", 3))),
                                    col = list(status = c("Control" =  "orange", "High" = "purple")))

png(file="Heatmap_mito.png", width=700, height=600, res=150)
print(Heatmap(z.mat, name = "z-score",
              col = myRamp,   
              show_row_name = TRUE,
              cluster_columns = FALSE,
              rect_gp = gpar(col = "darkgrey", lwd=0.5),
              top_annotation = ha1))
dev.off()
