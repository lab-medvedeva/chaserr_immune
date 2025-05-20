library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

my_pbmc_s <- readRDS("AIDA_data/ff5a921e-6e6c-49f6-9412-ad9682d23307.rds")
my_pbmc_s <- NormalizeData(my_pbmc_s)

# A. Umap with l1 annotation
Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$Annotation_Level1

p1 <- DimPlot(my_pbmc_s, group.by = "Annotation_Level1",
  label = TRUE, label.size = 2, repel = TRUE, raster=FALSE) + ggtitle(NULL) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        plot.tag = element_text(size = 8, face="plain")) +
    labs(title = "", tag = "A")

# B. Umap with CHASERR expression
Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$Annotation_Level2

p2 <- FeaturePlot(my_pbmc_s,
      features = c('ENSG00000272888'),
      label = FALSE,
      raster = FALSE) + scale_color_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        plot.tag = element_text(size = 8, face="plain")) +
  labs(title = "", tag = "B")

# C. Violin plot with CHASERR expression
Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$cell_type

p3 <- VlnPlot(
  my_pbmc_s,
  features = c('ENSG00000272888'),
  alpha = 0,
  sort = 'decreasing') + NoLegend() +
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
  labs(title = "", tag = "C")

g <- arrangeGrob(p1, p2, p3, ncol = 2, nrow = 7, layout_matrix= rbind(c(1,3), c(1,3), c(1,3), c(1,3), c(2,3), c(2,3), c(2,3)))
ggsave(paste0("chaserr_expression_main.png"), plot = g, units = "mm", width = 170, height = 110)
