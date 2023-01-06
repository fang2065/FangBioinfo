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

data <- read_tsv("Counts_TMT_Mito.txt", comment="#")
data<-data[!is.na(data$H1) & !is.na(data$C1),]

countdata <- data %>%
  as.data.frame() %>% 
  column_to_rownames("Accession") %>% # turn the geneid column into rownames
  as.matrix()

#box plot 
par(mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
boxplot(countdata[, 1:6],  ylim = c(0, 2.5), main="Boxplot normalized Intensities")

#limma test of significant protein 
#moderated P value calculated from https://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html
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


## Annoate the Accession id 
mouseUp <- UniProt.ws(10090) 
annot <- AnnotationDbi::select(x = mouseUp,
  columns = c("gene_primary"),
  keys = rownames(countdata),
  keytype = "UniProtKB")

annot <- annot[,-1] %>%
         rename(Accession=Entry, GeneID=Gene.Names..primary.) 

resannot <- as.data.frame(res.eb) %>% 
  rownames_to_column("Accession") %>% 
  left_join(annot, "Accession")

png(file="Volcanoplot_peptide_mito.png", width=800, height=800, res = 125)
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

resannotsort <- resannot[order(resannot$p.mod),]
write.table(resannotsort, file="Peptide_DE_mito.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

resannotsort_DN <- resannotsort %>% dplyr::filter(logFC <0.25 & p.mod < 0.01)
resannotsort_UP <- resannotsort %>% dplyr::filter(logFC >0.25 & p.mod < 0.01)

gene1_dn <- as.character(resannotsort_DN$GeneID)

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

ego_BP_dn_simplify_table <- simplify(ego_BP_dn, 
                                     cutoff=0.7, 
                                     by="p.adjust", 
                                     select_fun=min) %>%
                            as.data.frame()


resannotsort<-resannotsort %>% dplyr::filter(!is.na(GeneID) & !is.na(p.mod))
trend <-sapply(resannotsort$logFC, function(x){if(x>0) 'up' else 'down'})
resannot_order <- resannotsort[order(resannotsort$logFC, decreasing = T),]
pre_ranked_sig_genes<- data.frame(resannot_order, 'trend' = trend, 'rank' = 1:nrow(resannot_order), stringsAsFactors = F)
to_be_point_out <- rbind(pre_ranked_sig_genes_PV[pre_ranked_sig_genes_PV$Gene.symbol == "Cnp", ], 
                        pre_ranked_sig_genes_PV[pre_ranked_sig_genes_PV$Gene.symbol == "CD47", ])


GoID_mito = c('Tomm40',
            'Timm23',
            'Rhot2',
            'Rhot1',
            'Timm22')

to_be_point_out <- data.frame()
for(i in GoID_mito){
  tmp  <- pre_ranked_sig_genes[pre_ranked_sig_genes$GeneID == i,]
  to_be_point_out <- rbind(to_be_point_out, tmp)
}

tiff(file="fallplot_mito.tiff", width=1000, height=1000, res = 150)
ggplot(pre_ranked_sig_genes, aes(x=rank, y=logFC, color=logFC)) +
  geom_point(size=1)+
  geom_hline(yintercept = c(1.5,-1.5), linetype=2, size=0.25)+
  geom_hline(yintercept = c(0), linetype=1, size=0.5)+
  geom_vline(xintercept = 1636.5, linetype=2, size=0.25)+
  scale_color_gradient2(low="navy", high="firebrick3", mid="white", midpoint = 0)+
  geom_point(inherit.aes = F, data=to_be_point_out, aes(x=rank, y=logFC), size = 3, color = 'black')+
  geom_point(inherit.aes = F, data=to_be_point_out, aes(x=rank, y=logFC), size = 2, color = 'yellow')+
  ggrepel::geom_text_repel(inherit.aes = F, data = to_be_point_out, aes(x=rank, y=logFC, label=GeneID), size =5)+
  xlab('rank of differentially expressed proteins') +
  theme_bw()+
  theme(panel.grid = element_line(color = 'white'), legend.title.align = 0.5)
dev.off()

png(file="dotplot_peptide.png", width=2000, height=700, res = 120)
cowplot::plot_grid(Downregulated_pathway, Upregulated_pathway)
dev.off()


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


