library(Seurat)
library (ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(stringr)


# A. Cyclosporine scRNA-seq dataset clusters
counts_t_s <- readRDS("mouse2_T.rds")

p1 <- DimPlot(counts_t_s, reduction = "harmony_umap2", group.by = "harmony_clusters2", label = TRUE, label.size = 4, repel = TRUE, raster=FALSE) + 
  NoLegend() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "none",
        axis.title.y = element_text(size = 6),
        plot.title = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face="plain")) + 
  labs(title = "", tag = "A")

# C. Marker genes for annotation
p2 <- DotPlot(
  counts_t_s,
  features = c('Cd4','Cd8a','Foxp3','Gzmk','Igfbp4','Gzmb','Ctla2a','Mki67','Stmn1','Sostdc1','Cxcr5','Cd40lg','Ccr7','Cxcr6', 'Ifngr1','Il17a', 'Ccr2','Rorc'),
  group.by = "harmony_clusters2")+
  theme(text = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
  	axis.text.x = element_text(angle = 90, size = 6),
  	axis.text.y = element_text(size = 6),
  	plot.title = element_text(size = 6, hjust = 0.5),
  legend.key.size = unit(0.3, "cm"),  # Key size
  legend.text = element_text(size = 6),  # Text labels
  legend.title = element_text(size = 5), # Title
  plot.tag = element_text(size = 8, face="plain")) + 
  labs(title = "", tag = "B")

g <- arrangeGrob(p1, p2, ncol = 2, nrow = 1, layout_matrix= rbind(c(1), c(2)))
ggsave(paste0("mouse_cyclosporin_supplementary.png"), plot = g, units = "mm", width = 170, height = 150)

