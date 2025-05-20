library(Seurat)
library (ggplot2)
library(RColorBrewer)
library(gridExtra)

my_pbmc_s <- readRDS("AIDA_data/ff5a921e-6e6c-49f6-9412-ad9682d23307.rds") 
my_pbmc_s <- NormalizeData(my_pbmc_s)

# A. Umap with l1 annotation
Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$Annotation_Level1
p1 <- DimPlot(my_pbmc_s, group.by = "Annotation_Level1",
  label = TRUE, label.size = 2 , repel = TRUE, raster=FALSE) + NoLegend() + ggtitle(NULL) +
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
        plot.title = element_text(size = 6, face = "plain"),
        plot.tag = element_text(size = 8, face = "plain")) +
    labs(title = "", tag = "A")

# B. Umap with CHD2 expression
p2 <- FeaturePlot(my_pbmc_s,
      features = c('ENSG00000173575'),
      label = FALSE,
      raster = FALSE) + scale_color_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        legend.position="none",
        axis.title.x=element_blank(),#axis.text.x=element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        plot.tag = element_text(size = 8, face="plain")) +
  labs(title = "", tag = "B")

Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$cell_type

# C. Violin plot with CHD2 expression
p3 <- VlnPlot(
  my_pbmc_s,
  features = c('ENSG00000173575'),
  alpha = 0,
  sort = 'decreasing') + NoLegend() + coord_flip() +
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
        plot.tag = element_text(size = 8, face="plain")) + 
  labs(title = "", tag = "C")


# D. Wilcoxon test for T cell subtypes
data_name <- 'AIDA2'
marks <- read.csv(file = paste0(data_name, "_l3_wilcoxon_marks.csv"), sep="\t")
marks_adj <- marks
marks_adj$p_val_adj <-  signif(marks_adj$p_val_adj, digits=4)
chaserr_res <- dplyr::filter(marks_adj, marks_adj$gene %in%  c('ENSG00000272888', 'ENSG00000173575'))
chaserr_res[chaserr_res$gene == 'ENSG00000272888', 'gene'] = 'CHASERR'
chaserr_res[chaserr_res$gene == 'ENSG00000173575', 'gene'] = 'CHD2'
chaserr_res <- chaserr_res[order(-chaserr_res$avg_log2FC),]
chaserr_res$cluster= factor(chaserr_res$cluster, levels =  unique(chaserr_res$cluster))
p4 <- ggplot(chaserr_res, aes(x=gene , y = cluster, fill=avg_log2FC))+
theme_bw() +
geom_tile() +
scale_fill_distiller(palette = "RdBu", limits = c(-0.6, 0.6)) +
geom_text(data=chaserr_res, aes(x = gene, y = cluster, label = p_val_adj, angle = 0),  size=3) +
scale_x_discrete(drop = FALSE) +
scale_y_discrete(drop = FALSE) +
theme(panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      text = element_text(size = 6),
      axis.text.x = element_text(angle = 0, size = 6),
      axis.text.y = element_text(angle = 0, size = 6),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 6, face="plain"),
      plot.tag = element_text(size = 8, face="plain"),
      legend.position="bottom") +
labs(title = "", tag = "D")

g <- arrangeGrob(p1, p2, p3, p4, ncol = 2, nrow = 4, layout_matrix= rbind(c(1,3), c(2,3), c(4,4), c(4,4)))
ggsave(paste0("chd2_expression_supplementary.png"), plot = g, units = "mm", width = 170, height = 220)