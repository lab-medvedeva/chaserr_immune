library(Seurat)
library(BPCells)
library (ggplot2)
library(RColorBrewer)
library(gridExtra)

# CD4 T data https://doi.org/10.1016/j.xgen.2023.100473
# Convert to Seurat
data <- open_matrix_anndata_hdf5("res.cxg.h5ad")

 write_matrix_dir(
   mat = data,
   dir = "res.cxg"
 )

mat <- open_matrix_dir(dir = "res.cxg")

metadata <- read.csv("res.cxg.metadata", sep = '\t', row.names = 1) 
merged.object <- CreateSeuratObject(counts = mat, meta.data = metadata)
saveRDS(merged.object, "res.cxg.rds")

# DICE https://dice-database.org/
chaserr_dice <- read.table("AC013394.2_Expression_data.csv", header = F, quote = "\"")
celltype <- chaserr_dice[1, 'V1']
expr <- strsplit(chaserr_dice[1, 'V2'],",")[[1]]
chaserr_dice_melt <- data.frame(celltype, expr[2:length(expr)])
for(i in 2:nrow(chaserr_dice))
{
  celltype <- chaserr_dice[i, 'V1']
  expr <- strsplit(chaserr_dice[i, 'V2'],",")[[1]]
  df <- data.frame(celltype, expr[2:length(expr)])
  chaserr_dice_melt <- rbind(chaserr_dice_melt, df)
}
colnames(chaserr_dice_melt) <- c('celltype', 'value')
chaserr_dice_melt$gene <- 'CHASERR'
write.table(chaserr_dice_melt, file = paste0("chaserr_dice.csv"), sep="\t", col.names=TRUE, quote = FALSE)

chd2_dice <- read.table("CHD2_Expression_data.csv", header = F, quote = "\"")
celltype <- chd2_dice[1, 'V1']
expr <- strsplit(chd2_dice[1, 'V2'],",")[[1]]
chd2_dice_melt <- data.frame(celltype, expr[2:length(expr)])
for(i in 2:nrow(chd2_dice))
{
  celltype <- chd2_dice[i, 'V1']
  expr <- strsplit(chd2_dice[i, 'V2'],",")[[1]]
  df <- data.frame(celltype, expr[2:length(expr)])
  chd2_dice_melt <- rbind(chd2_dice_melt, df)
}
colnames(chd2_dice_melt) <- c('celltype', 'value')
chd2_dice_melt$gene <- 'CHD2'
write.table(chd2_dice_melt, file = paste0("chd2_dice.csv"), sep="\t", col.names=TRUE, quote = FALSE)

# A. CHASERR DICE
chaserr_dice_melt$celltype <- factor(chaserr_dice_melt$celltype, levels = unique(chaserr_dice_melt$celltype))
chaserr_dice_melt$value <- as.double(chaserr_dice_melt$value)

p1 <- ggplot(chaserr_dice_melt, aes(x = reorder(celltype, value, mean), y = value, fill = gene))+ 
  theme_bw() + 
  geom_boxplot(fatten = 0.8, outlier.shape = NA, alpha = 0.4) + coord_flip() +
  theme(
        strip.background = element_blank(),  
        strip.text.x = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        legend.title = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        plot.tag = element_text(size = 8, face="plain"),
        legend.position="none",
        axis.title = element_text(size = 6))+
  labs(title = "CHASERR", tag = "A")

# B. CHASERR CD4+ dataset L1
my_pbmc_s <- readRDS("res.cxg.rds")
my_pbmc_s <- NormalizeData(my_pbmc_s)
Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$cluster_L1
 
p2 <- VlnPlot(
  my_pbmc_s,
  features = c('LINC01578'),
  alpha = 0, sort = 'decreasing') + NoLegend() + 
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
        plot.tag = element_text(size = 8, face="plain")) + coord_flip()+
  labs(title = "", tag = "B")

# C. CHASERR CD4+ dataset L2
Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$cluster_L2
p3 <- VlnPlot(
  my_pbmc_s,
  features = c('LINC01578'),
  alpha = 0, sort = 'increasing') + NoLegend() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        plot.tag = element_text(size = 8, face="plain")) +
  labs(title = "", tag = "С")

# D. CHD2 DICE
chd2_dice_melt <- read.csv(paste0("chd2_dice.csv"), sep="\t")
chd2_dice_melt$celltype <- factor(chd2_dice_melt$celltype, levels = unique(chd2_dice_melt$celltype))
chd2_dice_melt$value <- as.double(chd2_dice_melt$value)

p4 <- ggplot(chd2_dice_melt, aes(x = reorder(celltype, value, mean), y = value, fill = gene))+ 
  theme_bw() + 
  geom_boxplot(fatten = 0.8, outlier.shape = NA, alpha = 0.4) + coord_flip() +
  theme(
        strip.background = element_blank(),  
        strip.text.x = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        legend.title = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        plot.tag = element_text(size = 8, face="plain"),
        legend.position="none",
        axis.title = element_text(size = 6)) +
  labs(title = "CHD2", tag = "D")

# CHD2 CD4+ dataset L1
Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$cluster_L1
 p5 <- VlnPlot(
  my_pbmc_s,
  features = c('CHD2'),
  alpha = 0, sort = 'decreasing') + NoLegend() + 
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
        plot.tag = element_text(size = 8, face="plain")) + coord_flip()+
  labs(title = "", tag = "E")

# CHD2 CD4+ dataset L2
Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$cluster_L2
p6 <- VlnPlot(
  my_pbmc_s,
  features = c('CHD2'),
  alpha = 0, sort = 'increasing') + NoLegend() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        plot.tag = element_text(size = 8, face="plain")) +
  labs(title = "", tag = "F")

g <- arrangeGrob(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 4, layout_matrix= rbind(c(1,2), c(3,3), c(4,5), c(6,6))) #p5,
ggsave(paste0("validation.png"), plot = g, units = "mm", width = 170, height = 220)


