library(Seurat)
library(BPCells)
library (ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(dplyr)

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

df <- read.csv('Tcells_B2M_CHASERR_CHD2_raw_Cq.csv', sep = '\t')
df_b2m <- df[df$target == 'b2m',]
df_chaserr <- df[df$target == 'chaserr',]
df_chaserr_naive <- df_chaserr[df_chaserr$sample == 'Naive','cq'] - df_b2m[df_b2m$sample == 'Naive', 'cq.mean'][1]
df_chaserr_cm <- df_chaserr[df_chaserr$sample == 'CM','cq'] - df_b2m[df_b2m$sample == 'CM', 'cq.mean'][1]
df_chaserr_th1 <- df_chaserr[df_chaserr$sample == 'Th1','cq'] - df_b2m[df_b2m$sample == 'Th1', 'cq.mean'][1]
df_chaserr_th17 <- df_chaserr[df_chaserr$sample == 'Th17','cq'] - df_b2m[df_b2m$sample == 'Th17', 'cq.mean'][1]
df_chaserr_th117 <- df_chaserr[df_chaserr$sample == 'Th1-17','cq'] - df_b2m[df_b2m$sample == 'Th1-17', 'cq.mean'][1]
df_chaserr_temra <- df_chaserr[df_chaserr$sample == 'TEMRA','cq'] - df_b2m[df_b2m$sample == 'TEMRA', 'cq.mean'][1]
df_chaserr_treg <- df_chaserr[df_chaserr$sample == 'Treg','cq'] - df_b2m[df_b2m$sample == 'Treg', 'cq.mean'][1]
df_chaserr_em <- df_chaserr[df_chaserr$sample == 'EM','cq'] - df_b2m[df_b2m$sample == 'EM', 'cq.mean'][1]
two_minus_ddct_chaserr_naive = 2**(-(df_chaserr_naive - df_chaserr_temra))
two_minus_ddct_chaserr_cm = 2**(-(df_chaserr_cm - df_chaserr_temra))
two_minus_ddct_chaserr_th1 = 2**(-(df_chaserr_th1 - df_chaserr_temra))
two_minus_ddct_chaserr_th17 = 2**(-(df_chaserr_th17 - df_chaserr_temra))
two_minus_ddct_chaserr_th117 = 2**(-(df_chaserr_th117 - df_chaserr_temra))
two_minus_ddct_chaserr_temra = 2**(-(df_chaserr_temra - df_chaserr_temra))
two_minus_ddct_chaserr_treg = 2**(-(df_chaserr_treg - df_chaserr_temra))
two_minus_ddct_chaserr_em = 2**(-(df_chaserr_em - df_chaserr_temra))
df_chd2 <- df[df$target == 'chd2',]
df_chd2_naive <- df_chd2[df_chd2$sample == 'Naive','cq'] - df_b2m[df_b2m$sample == 'Naive', 'cq.mean'][1]
df_chd2_cm <- df_chd2[df_chd2$sample == 'CM','cq'] - df_b2m[df_b2m$sample == 'CM', 'cq.mean'][1]
df_chd2_th1 <- df_chd2[df_chd2$sample == 'Th1','cq'] - df_b2m[df_b2m$sample == 'Th1', 'cq.mean'][1]
df_chd2_th17 <- df_chd2[df_chd2$sample == 'Th17','cq'] - df_b2m[df_b2m$sample == 'Th17', 'cq.mean'][1]
df_chd2_th117 <- df_chd2[df_chd2$sample == 'Th1-17','cq'] - df_b2m[df_b2m$sample == 'Th1-17', 'cq.mean'][1]
df_chd2_temra <- df_chd2[df_chd2$sample == 'TEMRA','cq'] - df_b2m[df_b2m$sample == 'TEMRA', 'cq.mean'][1]
df_chd2_treg <- df_chd2[df_chd2$sample == 'Treg','cq'] - df_b2m[df_b2m$sample == 'Treg', 'cq.mean'][1]
df_chd2_em <- df_chd2[df_chd2$sample == 'EM','cq'] - df_b2m[df_b2m$sample == 'EM', 'cq.mean'][1]

ddct_sd_chaserr = ((df_b2m$Cq.SD)**2 + (df_chaserr$Cq.SD)**2)**0.5
ddct_sd_chd2 = ((df_b2m$Cq.SD)**2 + (df_chd2$Cq.SD)**2)**0.5

two_minus_ddct_chd2_naive = 2**(-(df_chd2_naive - df_chd2_temra))
two_minus_ddct_chd2_cm = 2**(-(df_chd2_cm - df_chd2_temra))
two_minus_ddct_chd2_th1 = 2**(-(df_chd2_th1 - df_chd2_temra))
two_minus_ddct_chd2_th17 = 2**(-(df_chd2_th17 - df_chd2_temra))
two_minus_ddct_chd2_th117 = 2**(-(df_chd2_th117 - df_chd2_temra))
two_minus_ddct_chd2_temra = 2**(-(df_chd2_temra - df_chd2_temra))
two_minus_ddct_chd2_treg = 2**(-(df_chd2_treg - df_chd2_temra))
two_minus_ddct_chd2_em = 2**(-(df_chd2_em - df_chd2_temra))

cell_type <- c('Naive','Naive','Naive',
  'CM','CM','CM',
  'Th1','Th1','Th1',
  'Th17','Th17','Th17',
  'Th1-17','Th1-17','Th1-17',
  'TEMRA','TEMRA','TEMRA',
  'Treg','Treg','Treg',
  'EM','EM','EM')

relative_expression <- c(two_minus_ddct_chaserr_naive,
  two_minus_ddct_chaserr_cm,
  two_minus_ddct_chaserr_th1,
  two_minus_ddct_chaserr_th17,
  two_minus_ddct_chaserr_th117,
  two_minus_ddct_chaserr_temra,
  two_minus_ddct_chaserr_treg,
  two_minus_ddct_chaserr_em)
gene <- c('CHASERR',
  'CHASERR',
  'CHASERR',
  'CHASERR',
  'CHASERR',
  'CHASERR',
  'CHASERR',
  'CHASERR')
rrr1 <-  data.frame(cell_type, relative_expression, gene)

cell_type <- c('Naive','Naive','Naive',
  'CM','CM','CM',
  'Th1','Th1','Th1',
  'Th17','Th17','Th17',
  'Th1-17','Th1-17','Th1-17',
  'TEMRA','TEMRA','TEMRA',
  'Treg','Treg','Treg',
  'EM','EM','EM')
relative_expression <- c(two_minus_ddct_chd2_naive,
  two_minus_ddct_chd2_cm,
  two_minus_ddct_chd2_th1,
  two_minus_ddct_chd2_th17,
  two_minus_ddct_chd2_th117,
  two_minus_ddct_chd2_temra,
  two_minus_ddct_chd2_treg,
  two_minus_ddct_chd2_em)
gene <- c('CHD2',
  'CHD2',
  'CHD2',
  'CHD2',
  'CHD2',
  'CHD2',
  'CHD2',
  'CHD2')

rrr2 <-  data.frame(cell_type, relative_expression, gene)


#
df <- read.csv('admin_2025_08_28_16_46_25_CT030971_b2m_ppia_chasserr_chd2_immune.csv', sep = ';')
df <- df[order(df$Sample),]
df_chaserr <- df[df$Target == 'chaserr',]
df_chd2 <- df[df$Target == 'chd2',]
df_b2m <- df[df$Target == 'b2m',]
dct_chaserr = df_chaserr$Mean.Cq - df_b2m$Mean.Cq
dct_chd2 = df_chd2$Mean.Cq - df_b2m$Mean.Cq
dct_ppia = df_b2m$Mean.Cq - df_b2m$Mean.Cq
ddct_sd_chaserr = ((df_b2m$Cq.SD)**2 + (df_chaserr$Cq.SD)**2)**0.5
ddct_sd_chd2 = ((df_b2m$Cq.SD)**2 + (df_chd2$Cq.SD)**2)**0.5
two_minus_ddct_chaserr = 2**(-(dct_chaserr - dct_chaserr[df_chaserr$Sample == 'TEMRA']))
two_minus_ddct_chd2 = 2**(-(dct_chd2 - dct_chd2[df_chd2$Sample == 'TEMRA']))
df1 <- data.frame(two_minus_ddct_chaserr, ddct_sd_chaserr, df_chaserr$Sample, 'CHASERR')
colnames(df1) <- c('relative_expression', 'sd', 'cell_type', 'gene')
df2 <- data.frame(two_minus_ddct_chd2, ddct_sd_chd2, df_chd2$Sample, 'CHD2')
colnames(df2) <- c('relative_expression', 'sd', 'cell_type', 'gene')
df <- rbind(df1,df2)
symnum.args <- list(cutpoints = c(0, 0.00001, 0.0001, 0.001, 0.05, Inf), symbols = c("****", "***", "**", "*", ""))

comparisons = list(
c('Naive','EM'),
c('CM','EM'),
c('CM','TEMRA'),
c('Naive','TEMRA'))

p11 <- ggplot(df1, aes(x = cell_type, y = relative_expression, fill = cell_type)) +
  geom_hline(yintercept=1, color="red") +
  theme_bw() +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(x=cell_type, ymin=relative_expression-sd, ymax=relative_expression+sd), width=0.4, alpha=0.9) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x =  element_text(size = 6),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face="plain")) + 
    stat_compare_means(data = rrr1, method='t.test', label = "p.signif", comparisons = comparisons, symnum.args = symnum.args, size = 2) + # ref.group=".all.", 
  labs(title = "CHASERR", y = "Relative expression", tag = "A")

comparisons = list(c('Th1', 'Treg'))

p12 <- ggplot(df2, aes(x = cell_type, y = relative_expression, fill = cell_type)) +
  geom_hline(yintercept=1, color="red") +
  theme_bw() +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(x=cell_type, ymin=relative_expression-sd, ymax=relative_expression+sd), width=0.4, alpha=0.9) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x =  element_text(size = 6),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.2, "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face="plain")) + 
    stat_compare_means(data = rrr2, method='t.test', label = "p.signif", comparisons = comparisons, symnum.args = symnum.args, size = 2) + #ref.group=".all.", 
  labs(title = "CHD2", y = "Relative expression", tag = "")

p <- ggarrange(p11, p12, common.legend = TRUE, legend = 'right')

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

# CHASERR DICE
chaserr_dice_melt <- read.csv(paste0("chaserr_dice.csv"), sep="\t")
chaserr_dice_melt$celltype <- factor(chaserr_dice_melt$celltype, levels = unique(chaserr_dice_melt$celltype))
chaserr_dice_melt$value <- as.double(chaserr_dice_melt$value)

p2 <- ggplot(chaserr_dice_melt, aes(x = reorder(celltype, value, mean), y = value, fill = gene))+ 
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
  labs(title = "CHASERR", tag = "B")


# Differential expression results from DICE for CHASERR
chaserr_de <- read.csv('AC013394.2_DE_data.csv')
chaserr_copy <- chaserr_de
chaserr_copy$Cell1 <- chaserr_de$Cell2
chaserr_copy$Cell2 <- chaserr_de$Cell1
chaserr_copy$lfc <- -chaserr_copy$lfc
chaserr_de <- rbind(chaserr_de, chaserr_copy)
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

chaserr_de_wide <- dcast(chaserr_de, Cell1 ~ Cell2)
rownames(chaserr_de_wide) <- chaserr_de_wide$Cell1
chaserr_de_wide <- chaserr_de_wide[levels,c('Cell1', levels)]
chaserr_de_wide[upper.tri(chaserr_de_wide, diag=F)] <- NA 
chaserr_de <- melt(chaserr_de_wide)
colnames(chaserr_de) <- c('Cell1', 'Cell2', 'lfc')
chaserr_de$Cell1 <- factor(chaserr_de$Cell1, levels = rev(levels))
chaserr_de$Cell2 <- factor(chaserr_de$Cell2, levels = rev(levels))
p3 <- ggplot(chaserr_de, aes(x=Cell2, y = Cell1, fill=lfc))+
theme_bw() +
geom_tile() +
scale_fill_distiller(palette = "RdYlBu", limits = c(-2.3, 2.3), na.value = 'white') +# geom_text(data=chaserr_res, aes(x = names, y = group, label = p_val_adj_str, angle = 0), size=2.5) +
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
      legend.title=element_text(size = 6))  +
labs(title = "", tag = "C")

# CHASERR CD4+ dataset L2
my_pbmc_s <- readRDS("res.cxg.rds")
my_pbmc_s <- NormalizeData(my_pbmc_s)
Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$cluster_L2
p4 <- VlnPlot(
  my_pbmc_s,
  features = c('LINC01578'),
  alpha = 0,
  sort = 'decreasing') + NoLegend() + 
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
  labs(title = "", tag = "D")

g <- arrangeGrob(p, p2, p3, p4, ncol = 3, nrow = 7, layout_matrix= rbind(c(1,1,1), c(1,1,1), c(2,2,4),c(2,2,4), c(3,3,4), c(3,3,4), c(3,3,4)))
ggsave(paste0("validation.png"), plot = g, units = "mm", width = 170, height = 200)
