library(Seurat)
library (ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(stringr)

# A. CHASERR knockdown + cyclosporin 24h qPCR in human fibroblasts, 2nd experiment
df1 = read.csv("admin_2022-09-21 15-50-56_CT030971_CHASERR_KD_HOUSEKEEPING.csv", sep =',')
df1 = df1[c('Target', 'Sample', 'Cq')]

df2 = read.csv("admin_2022-09-21 19-05-24_CT030971_CHASERR_KD_CHASERR_CHD2.csv", sep =',')
df2 = df2[c('Target', 'Sample', 'Cq')]

df3 = read.csv("admin_2022-06-21 16-48-07_CT030971.csv", sep =',')
df3 = df3[c('Target', 'Sample', 'Cq')]

df4 = read.csv("admin_2022-06-21 14-56-12_CT030971.csv", sep =',')
df4 = df4[c('Target', 'Sample', 'Cq')]

expression = rbind(rbind(rbind(df1,df2),df3),df4)
expression$Sample <- str_split_fixed(expression$Sample, "-", 2)[,1]
cq_control_control= expression[expression$Sample == 362,]
cq_kd_control= expression[expression$Sample == 360,]
cq_control_cycl= expression[expression$Sample == 377,]
cq_kd_cycl= expression[expression$Sample == 375,]

dct_chaserr_control_control = mean(cq_control_control[cq_control_control['Target'] == 'CHASERR','Cq']) - mean(cq_control_control[cq_control_control['Target'] == 'PPIA','Cq'])
dct_chd2_control_control = mean(cq_control_control[cq_control_control['Target'] == 'CHD2','Cq']) - mean(cq_control_control[cq_control_control['Target'] == 'PPIA','Cq'])

dct_chaserr_kd_control = mean(cq_kd_control[cq_kd_control['Target'] == 'CHASERR','Cq']) - mean(cq_kd_control[cq_kd_control['Target'] == 'PPIA','Cq'])
dct_chd2_kd_control = mean(cq_kd_control[cq_kd_control['Target'] == 'CHD2','Cq']) - mean(cq_kd_control[cq_kd_control['Target'] == 'PPIA','Cq'])
two_minus_ddct_chaserr_kd_control = 2**(-(dct_chaserr_kd_control - dct_chaserr_control_control))
two_minus_ddct_chd2_kd_control = 2**(-(dct_chd2_kd_control - dct_chd2_control_control))

dct_chaserr_control_cycl = mean(cq_control_cycl[cq_control_cycl['Target'] == 'CHASERR','Cq']) - mean(cq_control_cycl[cq_control_cycl['Target'] == 'PPIA','Cq'])
dct_chd2_control_cycl = mean(cq_control_cycl[cq_control_cycl['Target'] == 'CHD2','Cq']) - mean(cq_control_cycl[cq_control_cycl['Target'] == 'PPIA','Cq'])
two_minus_ddct_chaserr_control_cycl = 2**(-(dct_chaserr_control_cycl - dct_chaserr_control_control))
two_minus_ddct_chd2_control_cycl = 2**(-(dct_chd2_control_cycl - dct_chd2_control_control))

dct_chaserr_kd_cycl = mean(cq_kd_cycl[cq_kd_cycl['Target'] == 'CHASERR','Cq']) - mean(cq_kd_cycl[cq_kd_cycl['Target'] == 'PPIA','Cq'])
dct_chd2_kd_cycl = mean(cq_kd_cycl[cq_kd_cycl['Target'] == 'CHD2','Cq']) - mean(cq_kd_cycl[cq_kd_cycl['Target'] == 'PPIA','Cq'])
two_minus_ddct_chaserr_kd_cycl = 2**(-(dct_chaserr_kd_cycl - dct_chaserr_control_control))
two_minus_ddct_chd2_kd_cycl = 2**(-(dct_chd2_kd_cycl - dct_chd2_control_control))


genes <- c('CHASERR','CHASERR','CHASERR','CHD2','CHD2','CHD2')
state <- c('KD_Control','KD_Cyclosporine','Control_Cyclosporine','KD_Control','KD_Cyclosporine','Control_Cyclosporine')
relative_expression <- c(two_minus_ddct_chaserr_kd_control, two_minus_ddct_chaserr_kd_cycl, two_minus_ddct_chaserr_control_cycl,
                         two_minus_ddct_chd2_kd_control, two_minus_ddct_chd2_kd_cycl, two_minus_ddct_chd2_control_cycl)
df <- data.frame(genes, state, relative_expression)
df$Sample <- paste0(df$genes, df$state)
df$Sample <- factor(df$Sample, levels = df$Sample)
p1 <- ggplot(df, aes(x=Sample, fill= state, y = relative_expression ))  +
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
         axis.ticks.x=element_blank(),
        legend.title = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        plot.title = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face="plain")) + 
  labs(title = "", tag = "A")

# B. Cyclosporine scRNA-seq dataset clusters
counts_t_s <- readRDS("mouse2_T.rds")

p2 <- DimPlot(counts_t_s, reduction = "harmony_umap2", group.by = "harmony_clusters2", label = TRUE, label.size = 4, repel = TRUE, raster=FALSE) + 
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
  labs(title = "", tag = "B")

# C. Marker genes for annotation
p3 <- DotPlot(
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
  labs(title = "", tag = "C")

g <- arrangeGrob(p1, p2, p3, ncol = 3, nrow = 3, layout_matrix= rbind(c(1,1,1), c(2,3,3), c(2,3,3)))
ggsave(paste0("mouse_cyclosporin_supplementary.png"), plot = g, units = "mm", width = 170, height = 150)

