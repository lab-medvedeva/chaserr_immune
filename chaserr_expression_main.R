library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggrepel)

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
        axis.title.y = element_blank(),
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
      pt.size = 0.01,
      raster = FALSE) + scale_color_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) +  #scale_colour_gradient2(low = "blue", mid = "yellow", high = "red") +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        plot.tag = element_text(size = 8, face="plain")) +
  labs(title = "", tag = "B")

# C. Wilcoxon test for T cell subtypes
data_name <- 'AIDA2'
marks <- read.csv(file = paste0(data_name, "_l3_wilcoxon_marks_pseudobulk.csv"), sep="\t")
marks_adj <- marks
chaserr_res <- dplyr::filter(marks_adj, marks_adj$gene %in%  c('ENSG00000272888')) #, 'ENSG00000173575'
chaserr_res[chaserr_res$gene == 'ENSG00000272888', 'gene'] = 'CHASERR'
chaserr_res[chaserr_res$gene == 'ENSG00000173575', 'gene'] = 'CHD2'
chaserr_res <- chaserr_res[order(-chaserr_res$avg_log2FC),]

chaserr_res <- dplyr::filter(chaserr_res, p_val_adj < 0.05) #, 'ENSG00000173575'

chaserr_res$p_val_adj <- -log(chaserr_res$p_val_adj, 10)

chaserr_res$significance <- NA
chaserr_res[chaserr_res$avg_log2FC < 0,'significance'] <- "Down"
chaserr_res[chaserr_res$avg_log2FC > 0,'significance'] <- "Up"
chaserr_res$significance <- factor(chaserr_res$significance, levels = c("Up", "Down"))

colors <- c( "Up" = "#E31A1C", "Down" =  "#1F78B4")

chaserr_res[chaserr_res$cluster == 'CD4+-T-naive', 'cluster'] = "CD4+ T naive"
chaserr_res[chaserr_res$cluster == 'CD8+-T-GZMBhi', 'cluster'] = "CD8+ T GZMBhi"
chaserr_res[chaserr_res$cluster == 'CD8+-T', 'cluster'] = "other CD8+ T"
chaserr_res[chaserr_res$cluster == 'MAIT', 'cluster'] = "MAIT"
chaserr_res[chaserr_res$cluster == 'CD8+-T-naive', 'cluster'] = "CD8+ T naive"
chaserr_res[chaserr_res$cluster == 'gdT-GZMBhi', 'cluster'] = "gdT GZMBhi"
chaserr_res[chaserr_res$cluster == 'CD4+-T-cyt', 'cluster'] = "CD4+ T cyt"
chaserr_res[chaserr_res$cluster == 'gdT-GZMKhi', 'cluster'] = "gdT GZMKhi"
chaserr_res[chaserr_res$cluster == 'T', 'cluster'] = "other T"
chaserr_res[chaserr_res$cluster == 'CD8+-T-GZMKhi', 'cluster'] = "CD8+ T GZMKhi"
chaserr_res[chaserr_res$cluster == 'CD4+-T', 'cluster'] = "other CD4+ T"
chaserr_res[chaserr_res$cluster == 'CD4+-T-em', 'cluster'] = "CD4+ Tem"
chaserr_res[chaserr_res$cluster == 'Treg', 'cluster'] = "Treg"
chaserr_res[chaserr_res$cluster == 'CD4+-T-cm', 'cluster'] = "CD4+ Tcm"

p3 <- ggplot(chaserr_res, aes(x = p_val_adj, y = avg_log2FC, color = significance, label = cluster)) +
  theme_bw() +
  geom_point(alpha = 1, size = 0.5) +
  geom_text_repel(hjust = 0, vjust = 0, size = 1.5) +
  xlim(0, 35) +
  ylim(-0.4, 0.4) +
  scale_color_manual(values = colors, name = "CHASERR") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +  
  theme(panel.border = element_blank(), #        panel.grid.major = element_blank(),        panel.grid.minor = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 5),
        legend.title=element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        plot.title = element_blank(),
        plot.tag = element_text(size = 8, face="plain")) +
  labs(title = "",
        x = "-log10(adjusted p-value)",
        y = "average log2FC",
        tag = "C")

# D. Violin plot with CHASERR expression
Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$cell_type

p4 <- VlnPlot(
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
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        plot.tag = element_text(size = 8, face="plain")) + coord_flip()+
  labs(title = "", tag = "D")


g <- arrangeGrob(p1, p2, p3, p4, ncol = 2, nrow = 11, layout_matrix= rbind(c(1,4), c(1,4), c(1,4), c(1,4), c(2,4), c(2,4),c(2,4), c(3,4), c(3,4), c(3,4), c(3,4)))
ggsave(paste0("chaserr_expression_main.png"), plot = g, units = "mm", width = 170, height = 160)





# D. Wilcoxon test for T cell subtypes
data_name <- 'AIDA2'
marks <- read.csv(file = paste0(data_name, "_l3_wilcoxon_marks.csv"), sep="\t")

marks_adj <- marks[marks$p_val_adj < 0.05,]
marks_adj$p_val_adj <- -log(marks_adj$p_val_adj, 10)

chaserr_res <- dplyr::filter(marks_adj, marks_adj$gene == 'ENSG00000272888')
chaserr_res[chaserr_res$gene == 'ENSG00000272888', 'gene'] = 'CHASERR'
chaserr_res <- chaserr_res[order(-chaserr_res$p_val_adj),]


chaserr_res[chaserr_res$cluster == 'CD4+_T_naive', 'cluster'] = "CD4+ T naive"
chaserr_res[chaserr_res$cluster == 'CD8+_T_GZMBhi', 'cluster'] = "CD8+ T GZMBhi"
chaserr_res[chaserr_res$cluster == 'CD8+_T', 'cluster'] = "other CD8+ T"
chaserr_res[chaserr_res$cluster == 'MAIT', 'cluster'] = "MAIT"
chaserr_res[chaserr_res$cluster == 'CD8+_T_naive', 'cluster'] = "CD8+ T naive"
chaserr_res[chaserr_res$cluster == 'gdT_GZMBhi', 'cluster'] = "gdT GZMBhi"
chaserr_res[chaserr_res$cluster == 'CD4+_T_cyt', 'cluster'] = "CD4+ T cyt"
chaserr_res[chaserr_res$cluster == 'gdT_GZMKhi', 'cluster'] = "gdT GZMKhi"
chaserr_res[chaserr_res$cluster == 'T', 'cluster'] = "other T"
chaserr_res[chaserr_res$cluster == 'CD8+_T_GZMKhi', 'cluster'] = "CD8+ T GZMKhi"
chaserr_res[chaserr_res$cluster == 'CD4+_T', 'cluster'] = "other CD4+ T"
chaserr_res[chaserr_res$cluster == 'CD4+_T_em', 'cluster'] = "CD4+ Tem"
chaserr_res[chaserr_res$cluster == 'Treg', 'cluster'] = "Treg"
chaserr_res$cluster= factor(chaserr_res$cluster, levels =  unique(chaserr_res$cluster))

p4 <- ggplot(chaserr_res, aes(x=cluster , y = p_val_adj, fill=avg_log2FC))+
theme_bw() +
geom_bar(position="dodge", stat="identity") +
scale_fill_distiller(palette = "RdBu", limits = c(-0.51, 0.51)) +
scale_x_discrete(drop = FALSE) +
scale_y_discrete(drop = FALSE) +
ylim(0, 280)+
theme(panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      text = element_text(size = 6),
      axis.text.x = element_text(angle = 20, size = 6),
      axis.text.y = element_text(angle = 0, size = 6),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(), 
      axis.title.y = element_text(angle = 90, size = 6),
      plot.title = element_text(size = 6, face="plain"),
      plot.tag = element_text(size = 8, face="plain"),
      legend.position="right") +
labs(title = "", y="-log10(p_val_adj)", tag = "D")
