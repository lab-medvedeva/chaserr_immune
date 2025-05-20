library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(clusterProfiler)
library("org.Mm.eg.db") 
 
# A. CHASERR knockdown + cyclosporin 20h qPCR in human fibroblasts, 1st experiment

df1 = read.csv("admin_2022-09-21 15-50-56_CT030971_CHASERR_KD_HOUSEKEEPING.csv", sep =',')
df1 = df1[c('Target', 'Sample', 'Cq')]

df2 = read.csv("admin_2022-09-21 19-05-24_CT030971_CHASERR_KD_CHASERR_CHD2.csv", sep =',')
df2 = df2[c('Target', 'Sample', 'Cq')]

df3 = read.csv("admin_2022-06-21 16-48-07_CT030971.csv", sep =',')
df3 = df3[c('Target', 'Sample', 'Cq')]

df4 = read.csv("admin_2022-06-21 14-56-12_CT030971.csv", sep =',')
df4 = df4[c('Target', 'Sample', 'Cq')]

expression = rbind(rbind(rbind(df1,df2),df3),df4)

cq_control_control= expression[expression$Sample == 343,]
cq_kd_control= expression[expression$Sample == 331,]
cq_control_cycl= expression[expression$Sample == 336,]
cq_kd_cycl= expression[expression$Sample == 324,]

dct_chaserr_control_control = cq_control_control[cq_control_control['Target'] == 'CHASERR','Cq'] - cq_control_control[cq_control_control['Target'] == 'PPIA','Cq']
dct_chd2_control_control = cq_control_control[cq_control_control['Target'] == 'CHD2','Cq'] - cq_control_control[cq_control_control['Target'] == 'PPIA','Cq']

dct_chaserr_kd_control = cq_kd_control[cq_kd_control['Target'] == 'CHASERR','Cq'] - cq_kd_control[cq_kd_control['Target'] == 'PPIA','Cq']
dct_chd2_kd_control = cq_kd_control[cq_kd_control['Target'] == 'CHD2','Cq'] - cq_kd_control[cq_kd_control['Target'] == 'PPIA','Cq']
two_minus_ddct_chaserr_kd_control = 2**(-(dct_chaserr_kd_control - dct_chaserr_control_control))
two_minus_ddct_chd2_kd_control = 2**(-(dct_chd2_kd_control - dct_chd2_control_control))

dct_chaserr_control_cycl = cq_control_cycl[cq_control_cycl['Target'] == 'CHASERR','Cq'] - cq_control_cycl[cq_control_cycl['Target'] == 'PPIA','Cq']
dct_chd2_control_cycl = cq_control_cycl[cq_control_cycl['Target'] == 'CHD2','Cq'] - cq_control_cycl[cq_control_cycl['Target'] == 'PPIA','Cq']
two_minus_ddct_chaserr_control_cycl = 2**(-(dct_chaserr_control_cycl - dct_chaserr_control_control))
two_minus_ddct_chd2_control_cycl = 2**(-(dct_chd2_control_cycl - dct_chd2_control_control))

dct_chaserr_kd_cycl = cq_kd_cycl[cq_kd_cycl['Target'] == 'CHASERR','Cq'] - cq_kd_cycl[cq_kd_cycl['Target'] == 'PPIA','Cq']
dct_chd2_kd_cycl = cq_kd_cycl[cq_kd_cycl['Target'] == 'CHD2','Cq'] - cq_kd_cycl[cq_kd_cycl['Target'] == 'PPIA','Cq']
two_minus_ddct_chaserr_kd_cycl = 2**(-(dct_chaserr_kd_cycl - dct_chaserr_control_control))
two_minus_ddct_chd2_kd_cycl = 2**(-(dct_chd2_kd_cycl - dct_chd2_control_control))

genes <- c('CHASERR','CHASERR','CHASERR','CHD2','CHD2','CHD2')
state <- c('KD_Control','KD_Cyclosporine','Control_Cyclosporine','KD_Control','KD_Cyclosporine','Control_Cyclosporine')

relative_expression <- c(two_minus_ddct_chaserr_kd_control, two_minus_ddct_chaserr_kd_cycl, two_minus_ddct_chaserr_control_cycl,
                         two_minus_ddct_chd2_kd_control, two_minus_ddct_chd2_kd_cycl, two_minus_ddct_chd2_control_cycl)

df <- data.frame(genes, state, relative_expression)
df$Sample <- paste0(df$genes, df$state)
df$Sample <- factor(df$Sample, levels = df$Sample)
p0 <- ggplot(df, aes(x=Sample, fill= state, y = relative_expression ))  +
  facet_grid(~ genes, scales = "free", space = "free" ) + 
  theme_bw() +
  geom_hline(yintercept=1, color="red") +
  geom_col(position = 'dodge') +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face="plain")) + 
  labs(title = "", tag = "A")

# B. Cyclosporine scRNA-seq dataset UMAP knockdown
counts_t_s <- readRDS("mouse2_T.rds")
counts_t_s@meta.data$cells <- rownames(counts_t_s@meta.data)
counts_t_s <- subset(counts_t_s, cells %in% counts_t_s@meta.data[!is.na(counts_t_s@meta.data$celltype_l2),'cells'] )

group_cell <- counts_t_s@meta.data$sample
group_cell[group_cell %in% c('EAU1','EAU2','EAU3','EAU4')] <- 'EAU'
group_cell[group_cell %in% c('CTL1','CTL2','CTL3','CTL4')] <- 'CTL'
group_cell[group_cell %in% c('CSA1','CSA2','CSA3','CSA4')] <- 'CSA'

counts_t_s@meta.data$de_group <- group_cell

p1 <- DimPlot(counts_t_s, reduction = "harmony_umap2", group.by = "celltype_l2", label = TRUE,
 label.size = 2, repel = TRUE, raster=FALSE, order = c("CD8 CTL", "CD4 other", "Pro T", "CD4 Th1", "CD4 Tfh", "CD4 Th17", "CD8 naive", "Treg", "CD4 naive")) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        plot.title = element_text(size = 6, face = "plain"),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face = "plain")) + 
labs(title = "", tag = "B")

# C. Cyclosporine scRNA-seq dataset Chaserr expression
counts_t_s@meta.data$cells <- rownames(counts_t_s@meta.data)
counts_t_s <- subset(x = counts_t_s, cells %in% counts_t_s@meta.data[!is.na(counts_t_s@meta.data$celltype_l2),'cells'])
subset_cells <- WhichCells(counts_t_s, expression = Chaserr > 0, accept.low = 0.1, use.raw=T)
counts_t_s_filtered <- subset(x = counts_t_s, cells %in% subset_cells)

p2 <- FeaturePlot(counts_t_s_filtered,
      features = c('Chaserr'),
      reduction = "harmony_umap2",
      label = FALSE,
      raster = FALSE) + scale_color_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        plot.title = element_text(size = 6, face = "plain"),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face = "plain")) + 
labs(title = "", tag = "C")

# D. Chaserr expression violin plot
p3 <- VlnPlot(
  counts_t_s_filtered,
  features = c('Chaserr'),
  group.by = 'celltype_l2',
  alpha = 0, sort = 'increasing') + 
  stat_summary(fun = median, geom='point', size = 10, colour = "black", shape = 95) +
  NoLegend() +
  #split.by = "de_group"
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face="plain")) + 
  labs(title = "", tag = "D")

# E. GO:BP for genes with positive correlation
res_adj1 <- read.csv(file = paste0("mouse_cyclosporin_corr.csv"), sep="\t")
res_adj1_pos <- res_adj1[res_adj1$corrs > 0,]
ego_BP_pos <- enrichGO(gene      = res_adj1_pos$genes,
                   OrgDb         = org.Mm.eg.db,
                   universe      = rownames(counts_t_s_filtered),
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05)

p4 <- cnetplot(ego_BP_pos, showCategory = 5, categorySize="pvalue", foldChange=NULL, cex_label_category = 0.4, cex_label_gene = 0.4) + 
labs(title = paste0("Chaserr Spearman's rho > 0"), tag = "E") + 
  theme(axis.title.x=element_blank(), 
        text = element_text(size = 6), 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 6, face="plain"),
        plot.tag = element_text(size = 8, face="plain"),
        legend.position = 'none')

# F.,G. Wilcoxon test results for Chaserr and Chd2
marks <- read.csv(paste0("marks_CTL_EAU_CSA.csv"), sep="\t")
chaserr_res <- dplyr::filter(marks, gene == 'Chaserr')
chaserr_res$comp= factor(chaserr_res$comp, levels = c('CTL_EAU', 'EAU_CSA'))
chaserr_res$celltype= factor(chaserr_res$celltype, levels =  unique(chaserr_res$celltype))
chaserr_res <- chaserr_res[order(-chaserr_res$avg_log2FC),]
p6 <- ggplot(chaserr_res, aes(x = celltype, y = comp , fill=avg_log2FC)) +
  theme_bw() +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", limits = c(-1.5,1.5)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face="plain")) +
  labs(title = "Chaserr", tag = "F")

Chd2_res <- dplyr::filter(marks, gene == 'Chd2')
Chd2_res$comp= factor(Chd2_res$comp, levels = c('CTL_EAU', 'EAU_CSA'))
Chd2_res$celltype= factor(Chd2_res$celltype, levels =  unique(Chd2_res$celltype))
Chd2_res <- Chd2_res[order(-Chd2_res$avg_log2FC),]
p7 <- ggplot(Chd2_res, aes(x = celltype, y = comp, fill=avg_log2FC))+
theme_bw() +
geom_tile() +
scale_fill_distiller(palette = "RdBu", limits = c(-1.4,1.4)) +
scale_x_discrete(drop = FALSE) +
scale_y_discrete(drop = FALSE) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face="plain")) +
labs(title = "Chd2", tag = "G")

g <- arrangeGrob(p0, p1, p2, p3, p4, p6, p7, ncol = 2, nrow = 5, layout_matrix= rbind(c(1,1),c(2,3), c(4,4), c(5,5), c(6,7))) 
ggsave(paste0("mouse_cyclosporin.png"), plot = g, units = "mm", width = 170, height = 225)
