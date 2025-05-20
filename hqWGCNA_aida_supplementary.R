library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(hdWGCNA)

my_pbmc_s <- readRDS('aida_metacells.rds')
obj <- my_pbmc_s@misc$tutorial$wgcna_metacell_obj
obj[["RNA"]]$data <- obj[["RNA"]]$counts
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:30)
obj <- subset(x = obj, subset = cell_type !=  'CD4-positive, alpha-beta T cell')
obj <- subset(x = obj, subset = cell_type !=  'CD8-positive, alpha-beta T cell')
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

Idents(object = obj) <- obj@meta.data$cell_type

p1 <- VlnPlot(
  obj,
  features = c('ENSG00000272888'),
  alpha = 0, sort = 'decreasing') + NoLegend() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        plot.tag = element_text(size = 8, face="plain")) + coord_flip()+
  labs(title = "CHASERR", tag = "")

p2 <- VlnPlot(
  obj,
  features = c('ENSG00000173575'),
  alpha = 0, sort = 'decreasing') + NoLegend() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        plot.tag = element_text(size = 8, face="plain")) + coord_flip()+
  labs(title = "CHD2", tag = "")

g <- arrangeGrob(p1, p2, ncol = 2, nrow = 1, layout_matrix= rbind(c(1,2)))
ggsave(paste0("metacells_supplementary.png"), plot = g, units = "mm", width = 170, height = 70)