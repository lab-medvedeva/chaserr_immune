library(Seurat)
library(BPCells)
library (ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(dplyr)

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

# A. CHD2 DICE
chd2_dice_melt <- read.csv(paste0("chd2_dice.csv"), sep="\t")
chd2_dice_melt$celltype <- factor(chd2_dice_melt$celltype, levels = unique(chd2_dice_melt$celltype))
chd2_dice_melt$value <- as.double(chd2_dice_melt$value)

p1 <- ggplot(chd2_dice_melt, aes(x = reorder(celltype, value, mean), y = value, fill = gene))+ 
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
  labs(title = "CHD2", tag = "A")

# B. Differential expression results from DICE for CHD2
chd2_de <- read.csv('CHD2_DE_data.csv')
chd2_copy <- chd2_de
chd2_copy$Cell1 <- chd2_de$Cell2
chd2_copy$Cell2 <- chd2_de$Cell1
chd2_copy$lfc <- -chd2_copy$lfc
chd2_de <- rbind(chd2_de, chd2_copy)
levels = c('T cell, CD4, naive TREG',
'B cell, naive',
'T cell, CD4, naive',
'T cell, CD4, TH1/17',
'T cell, CD4, TH2',
'T cell, CD8, naive',
'T cell, CD4, TH17',
'Monocyte, non-classical',
'T cell, CD4, TH1',
'T cell, CD4, TFH',
'T cell, CD4, memory TREG',
'Monocyte, classical',
'NK cell, CD56dim CD16+',
'T cell, CD4, naive [activated]',
'T cell, CD8, naive [activated]')

chd2_de_wide <- dcast(chd2_de, Cell1 ~ Cell2)
rownames(chd2_de_wide) <- chd2_de_wide$Cell1
chd2_de_wide <- chd2_de_wide[levels,c('Cell1', levels)]
chd2_de_wide[upper.tri(chd2_de_wide, diag=F)] <- NA 
chd2_de <- melt(chd2_de_wide)
colnames(chd2_de) <- c('Cell1', 'Cell2', 'lfc')

chd2_de$Cell1 <- factor(chd2_de$Cell1, levels = rev(levels))
chd2_de$Cell2 <- factor(chd2_de$Cell2, levels = rev(levels))
p2 <- ggplot(chd2_de, aes(x=Cell2, y = Cell1, fill=lfc))+
theme_bw() +
geom_tile() +
scale_fill_distiller(palette = "RdYlBu", limits = c(-1.6, 1.6), na.value = 'white') +# geom_text(data=chaserr_res, aes(x = names, y = group, label = p_val_adj_str, angle = 0), size=2.5) +
scale_x_discrete(drop = FALSE) +
scale_y_discrete(drop = FALSE) +
theme(panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      text = element_text(size = 6),
      axis.text.x = element_text(angle = 90, size = 6),
      axis.text.y = element_text(angle = 0, size = 6),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 6, face="plain"),
      plot.tag = element_text(size = 8, face="plain"),
      legend.position="right",
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5),
      legend.title=element_text(size = 6)) +
labs(title = "", tag = "B")

# C. CHD2 CD4+ dataset L2
my_pbmc_s <- readRDS("res.cxg.rds")
my_pbmc_s <- NormalizeData(my_pbmc_s)
Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$cluster_L2
p3 <- VlnPlot(
  my_pbmc_s,
  features = c('CHD2'),
  alpha = 0, sort = 'decreasing') + NoLegend() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        plot.tag = element_text(size = 8, face="plain")) + coord_flip()+
  labs(title = "", tag = "C")

# D. Wilcoxon test for T cell subtypes
data_name <- 'Validation'
marks <- read.csv(file = paste0(data_name, "_l2_wilcoxon_marks.csv"), sep="\t")
marks_adj <- marks
chaserr_res <- dplyr::filter(marks_adj, marks_adj$names %in%  c('LINC01578', 'CHD2'))
chaserr_res[chaserr_res$names == 'LINC01578', 'names'] = 'CHASERR'
chaserr_res[chaserr_res$names == 'CHD2', 'names'] = 'CHD2'
chaserr_res <- chaserr_res[order(-chaserr_res$logfoldchanges),]

chaserr_res$p_val_adj_str <- NA
chaserr_res$p_val_adj_str <- as.character(signif(chaserr_res$pvals_adj, digits=4))
chaserr_res[chaserr_res$pvals_adj == 0, 'p_val_adj_str'] = " < 2.225074e-308"
chaserr_res$group= factor(chaserr_res$group, levels =  unique(chaserr_res$group))

p4 <- ggplot(chaserr_res, aes(x=names, y = group, fill=logfoldchanges))+
theme_bw() +
geom_tile() +
scale_fill_distiller(palette = "RdBu", limits = c(-0.6, 0.6)) +
geom_text(data=chaserr_res, aes(x = names, y = group, label = p_val_adj_str, angle = 0), size=2.5) +
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
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 6, face="plain"),
      plot.tag = element_text(size = 8, face="plain"),
      legend.position="right",
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5),
      legend.title=element_text(size = 6)) +
labs(title = "", tag = "D")

g <- arrangeGrob(p1, p2, p3, p4, ncol = 3, nrow = 7, layout_matrix= rbind(c(1,1,3),c(1,1,3), c(2,2,3), c(2,2,3), c(2,2,3), c(4,4,4), c(4,4,4)))
ggsave(paste0("validation_supplementary.png"), plot = g, units = "mm", width = 170, height = 200)

