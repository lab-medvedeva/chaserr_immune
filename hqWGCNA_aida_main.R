library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(clusterProfiler)
library("org.Hs.eg.db")
library(BiocParallel)
library(enrichplot)
library(ggh4x)
library(ggpubr)

library(GO.db)
GO <- as.list(GOTERM)

enableWGCNAThreads(nThreads = 16)

my_pbmc_s <- readRDS('aida_metacells.rds')

obj <- my_pbmc_s@misc$tutorial$wgcna_metacell_obj
obj <- subset(x = obj, subset = cell_type !=  'CD4-positive, alpha-beta T cell')
obj <- subset(x = obj, subset = cell_type !=  'CD8-positive, alpha-beta T cell')
obj[["RNA"]]$data <- obj[["RNA"]]$counts
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:30)

obj@meta.data$cell_type[obj@meta.data$cell_type == "CD4-positive, alpha-beta cytotoxic T cell"] = 'CD4+ T cytotoxic'
obj@meta.data$cell_type[obj@meta.data$cell_type == "CD4-positive, alpha-beta T cell"] = 'CD4+ T other'
obj@meta.data$cell_type[obj@meta.data$cell_type == "CD8-positive, alpha-beta cytotoxic T cell"] = 'CD8+ T cytotoxic'
obj@meta.data$cell_type[obj@meta.data$cell_type == "CD8-positive, alpha-beta memory T cell"] = 'CD8+ T memory'
obj@meta.data$cell_type[obj@meta.data$cell_type == "CD8-positive, alpha-beta T cell"] = 'CD8+ T other'
obj@meta.data$cell_type[obj@meta.data$cell_type == "central memory CD4-positive, alpha-beta T cell"] = 'CD4+ TCM'
obj@meta.data$cell_type[obj@meta.data$cell_type == "effector memory CD4-positive, alpha-beta T cell"] = 'CD4+ TEM'
obj@meta.data$cell_type[obj@meta.data$cell_type == "naive thymus-derived CD4-positive, alpha-beta T cell"] = 'CD4+ T naive'
obj@meta.data$cell_type[obj@meta.data$cell_type == "naive thymus-derived CD8-positive, alpha-beta T cell"] = 'CD8+ T naive'
obj@meta.data$cell_type[obj@meta.data$cell_type == "regulatory T cell"] = 'CD4+ Treg'

# metacell UMAP
p1 <- DimPlot(obj, group.by = "cell_type",
  label = TRUE, label.size = 2, repel = TRUE, raster=FALSE) +  ggtitle(NULL) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        plot.title = element_text(size = 6, face = "plain"),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face = "plain"),
        legend.position = "none") +
  guides(color=guide_legend(ncol = 1)) +
  labs(title = "", tag = "A")

# Chaserr expression in metacells
p2 <- FeaturePlot(obj,
      features = c('ENSG00000272888'),
      label = FALSE,
      raster = FALSE) + scale_color_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) +
  theme(text = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    plot.title = element_text(size = 6, face = "plain"),
    axis.title = element_text(size = 6),
    plot.tag = element_text(size = 8, face = "plain")) +
labs(title = "", tag = "B")

# GO:BP plot
ego_BP_neg <- readRDS(paste0("effector memory CD4-positive, alpha-beta T cell_go_bp_neg.rds"))
p3 <-  barplot(ego_BP_neg, showCategory = 5 ) + labs(title = paste0("CD4+ TEM \n Spearman's rho < 0"), tag = "C") + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(), 
        text = element_text(size = 6), 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face="plain")) + 
  scale_fill_gradientn(colours=c('#FD7446', '#FFFFFF', '#709AE1'), limits=c(0, 0.05))

# FANTOM6 CHASERR knockdown signature scores
aso_df <- read.csv("DESeq2_genes_ASO_all.tsv", sep = '\t')
aso_df <-  dplyr::filter(aso_df, geneID %in% rownames(obj))
aso_df_fdr <- dplyr::filter(aso_df, fdr < 0.05)
aso_df_07 <- dplyr::filter(aso_df_fdr, perturb_id == 'ASO_G0272888_AD_07')
aso_df_07 <- aso_df_07[aso_df_07$geneID != 'ENSG00000272888',]
aso_df_10 <- dplyr::filter(aso_df_fdr, perturb_id == 'ASO_G0272888_AD_10')
aso_df_10 <- aso_df_10[aso_df_10$geneID != 'ENSG00000272888',]

markers <- list()
markers$aso_df_07_up <- unique( dplyr::filter(aso_df_07, log2FC >= 0 )$geneID)
markers$aso_df_07_down <- unique( dplyr::filter(aso_df_07, log2FC < 0 )$geneID)
markers$aso_df_10_up <- unique( dplyr::filter(aso_df_10, log2FC >= 0 )$geneID)
markers$aso_df_10_down <- unique( dplyr::filter(aso_df_10, log2FC < 0 )$geneID)

obj <- AddModuleScore(obj, features = markers, name = "add_module_score")

p11 <- VlnPlot(obj, features = paste0("add_module_score",1), group.by = "cell_type", alpha = 0, raster=FALSE) +
  theme(text = element_text(size = 8),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = 8),
  plot.title = element_text(size = 6, face="plain"),
  axis.title = element_text(size = 6),
  plot.tag = element_text(size = 8, face="plain")) +
  labs(title = "ASO_G0272888_AD_07, log2FC \u2265 0", tag = "E")
p12 <- VlnPlot(obj, features = paste0("add_module_score",2), group.by = "cell_type", alpha = 0, raster=FALSE) +
  theme(text = element_text(size = 8),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = 8),
  plot.title = element_text(size = 6, face="plain"),
  axis.title = element_text(size = 6),
  plot.tag = element_text(size = 8, face="plain")) +
  labs(title = "ASO_G0272888_AD_07, log2FC < 0", tag = "")
p13 <- VlnPlot(obj, features = paste0("add_module_score",3), group.by = "cell_type", alpha = 0, raster=FALSE) +
  theme(text = element_text(size = 8),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = 8),
  plot.title = element_text(size = 6, face="plain"),
  axis.title = element_text(size = 6),
  plot.tag = element_text(size = 8, face="plain")) +
  labs(title = "ASO_G0272888_AD_10, log2FC \u2265 0", tag = "")
p14 <- VlnPlot(obj, features = paste0("add_module_score",4), group.by = "cell_type", alpha = 0, raster=FALSE) +
  theme(text = element_text(size = 8),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = 8),
  plot.title = element_text(size = 6, face="plain"),
  axis.title = element_text(size = 6),
  plot.tag = element_text(size = 8, face="plain")) +
  labs(title = "ASO_G0272888_AD_10, log2FC < 0", tag = "")

p5 <- ggarrange(p11, p12, p13, p14, common.legend = TRUE, legend = 'bottom', nrow = 2, ncol = 2)

### GSVA umap
go_res <- readRDS('gsva_res_all.rds')
obj[["GSVA"]] <- CreateAssayObject(counts = go_res)
DefaultAssay(obj) <- "GSVA"
pattern_name = 'GO:0042110'
anno <- AnnotationDbi::select(GO.db, keys=c(pattern_name), keytype="GOID", columns=c("TERM","ONTOLOGY") )
res_adj1 <- read.csv("gsva_all_corr.csv", sep="\t")

p4 <-FeaturePlot(obj, pattern_name, pt.size = 0.5, reduction = "umap") + scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")))  +
 theme(axis.title.x=element_blank(),#axis.text.x=element_blank(),
        text = element_text(size = 6), #, face = "bold"
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face="plain")) +
labs(title = paste0(anno[anno$GOID == pattern_name, 'TERM'], "\n corr = ", res_adj1[res_adj1$genes == pattern_name,'corrs']), tag = "D")

g <- arrangeGrob(p1, p2, p3, p4, p5, ncol = 2, nrow = 4, layout_matrix= rbind(c(1,2), c(3,4), c(5,5))) 
ggsave(paste0("metacells.png"), plot = g,  units = "mm", width = 170, height = 225)

