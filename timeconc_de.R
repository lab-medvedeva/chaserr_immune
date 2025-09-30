library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggrepel)
library(ggpubr)

df1 <- read.csv('GSE90569_raw_counts_GRCh38.p13_NCBI.tsv', sep = '\t', row.names = 1)

treg_metadata <- read.csv('../treg_metadata.csv', sep = '\t')
treg_metadata <-  dplyr::filter(treg_metadata, time %in% c(0,1,2) )

samples1 <- treg_metadata
samples1$time <- as.character(samples1$time)
samples1$time <- factor(samples1$time, levels = c("0", "1", "2"))

df1 <- df1[,samples1$sample]
dds <- DESeqDataSetFromMatrix(countData = df1, design= ~ time , colData = samples1)
dds <- DESeq(dds)

res1 <- results(dds, name=, contrast=c("time", "1", "0"))
res1 <- res1[c('1106', '100507217'),]
res1$gene <- rownames(res1)
res1$time <- "0h_vs_1h"

res2 <- results(dds, name=, contrast=c("time", "2", "0"))
res2 <- res2[c('1106', '100507217'),]
res2$gene <- rownames(res2)
res2$time <- "0h_vs_2h"

res_GSE90569 <- rbind(res1,res2)
res_GSE90569$dataset <- "GSE90569"

df2 <- read.csv('GSE96538_raw_counts_GRCh38.p13_NCBI.tsv', sep = '\t', row.names = 1)

treg_metadata <- read.csv('../treg2_metadata.csv', sep = '\t')
treg_metadata <-  dplyr::filter(treg_metadata, time %in% c(0,1,2) )

samples2 <- treg_metadata
samples2$time <- as.character(samples2$time)
samples2$time <- factor(samples2$time, levels = c("0", "1", "2"))



df2 <- df2[,samples2$sample]
dds <- DESeqDataSetFromMatrix(countData = df2, design= ~ time , colData = samples2)
dds <- DESeq(dds)

res1 <- results(dds, contrast=c("time", "1", "0"))
res1 <- res1[c('1106', '100507217'),]
res1$gene <- rownames(res1)
res1$time <- "0h_vs_1h"

res2 <- results(dds, contrast=c("time", "2", "0"))
res2 <- res2[c('1106', '100507217'),]
res2$gene <- rownames(res2)
res2$time <- "0h_vs_2h"

res_GSE96538 <- rbind(res1,res2)
res_GSE96538$dataset <- "GSE96538"




df3 <- read.csv('GSE94396_raw_counts_GRCh38.p13_NCBI.tsv', sep = '\t', row.names = 1)
samples3 <- read.csv('../GSE94396_series_matrix.txt', sep = '\t')

colnames(samples3) <- c('sample', 'cell_type', 'time')
samples3 <-  dplyr::filter(samples3, time %in% c(0,2) )

samples3$time <- as.character(samples3$time)
samples3$time <- factor(samples3$time, levels = c("0", "2"))

df3 <- df3[,samples3$sample]
dds <- DESeqDataSetFromMatrix(countData = df3, design= ~ time , colData = samples3)
dds <- DESeq(dds)
res <- results(dds, name=, contrast=c("time", "2", "0"))

res <- res[c('1106', '100507217'),]
res$gene <- rownames(res)
res$time <- "0h_vs_2h"

res_GSE94396 <- res
res_GSE94396$dataset <- "GSE94396"


df4 <- read.csv('GSE52260_raw_counts_GRCh38.p13_NCBI.tsv', sep = '\t', row.names = 1)
samples4 <- read.csv('../th17_metadata_human.csv', sep = '\t')

samples4[samples4$time == '0h', 'time'] = '0'
samples4[samples4$time == '1h', 'time'] = '1'
samples4[samples4$time == '2h', 'time'] = '2'
samples4 <-  dplyr::filter(samples4, time %in% c('0','1','2') )
samples4$time <- factor(samples4$time, levels = c("0", "1", "2"))
df4 <- df4[,samples4$sample]
dds <- DESeqDataSetFromMatrix(countData = df4, design= ~ time , colData = samples4)
dds <- DESeq(dds)

res1 <- results(dds, contrast=c("time", "1", "0"))
res1 <- res1[c('1106', '100507217'),]
res1$gene <- rownames(res1)
res1$time <- "0h_vs_1h"

res2 <- results(dds, contrast=c("time", "2", "0"))
res2 <- res2[c('1106', '100507217'),]
res2$gene <- rownames(res2)
res2$time <- "0h_vs_2h"

res_GSE52260 <- rbind(res1,res2)
res_GSE52260$dataset <- "GSE52260"



df5 <- read.csv('../GSE140244_rnaseq_gene_counts.txt', sep = '\t', row.names = 1)
samples5 <- read.csv('../GSE140244_rnaseq_meta_data.txt', sep = '\t')
samples5 <- samples5[,c('ExpressionMatrix_SampleID', 'Time_point')]

colnames(samples5)<- c('sample', 'time')

samples5$time <- as.character(samples5$time)
samples5 <-  dplyr::filter(samples5, time %in% c('0', '2') )
samples5$time <- factor(samples5$time, levels = c("0", "2"))


df5 <- df5[,samples5$sample]
dds <- DESeqDataSetFromMatrix(countData = df5, design= ~ time , colData = samples5)
dds <- DESeq(dds)

res2 <- results(dds, contrast=c("time", "2", "0"))
res2 <- res2[c('ENSG00000173575', 'ENSG00000272888'),]
res2$gene <- rownames(res2)
res2$time <- "0h_vs_2h"

res_GSE140244 <- res2
res_GSE140244$dataset <- "GSE140244"



res <- rbind(res_GSE90569, res_GSE96538, res_GSE94396, res_GSE52260, res_GSE140244)
res$padj <- as.character(signif(res$padj, digits=4))



res[res$gene == '1106', 'gene'] = 'CHD2'
res[res$gene == '100507217', 'gene'] = 'CHASERR'

res[res$gene == 'ENSG00000173575', 'gene'] = 'CHD2'
res[res$gene == 'ENSG00000272888', 'gene'] = 'CHASERR'

res$gene <- factor(res$gene, levels = c('CHD2', 'CHASERR'))
res1 <- res[res$time == "0h_vs_1h",]
res1 <- res1[res1$dataset %in% c("GSE90569", "GSE96538", "GSE52260"),]
res1$dataset = factor(res1$dataset, levels = c("GSE90569", "GSE96538", "GSE52260"))

p1 <- ggplot(res1, aes(x=dataset , y = gene , fill=log2FoldChange))+
theme_bw() +
geom_tile() +
scale_fill_distiller(palette = "RdBu", limits = c(-1.60, 1.60)) +
geom_text(data=res1, aes(x = dataset, y = gene, label = padj, angle = 0), size=2.0) +
scale_x_discrete(drop = FALSE) +
scale_y_discrete(drop = FALSE) +
theme(panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      text = element_text(size = 5),
      axis.text.x = element_text(angle = 0, size = 5),
      axis.text.y = element_text(angle = 0, size = 5),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 6, face="plain"),
      plot.tag = element_text(size = 8, face="plain"),
      legend.position="right",
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5),
      legend.title=element_text(size = 5)) +
labs(title = "0h vs 1h", tag = "A")


res2 <- res[res$time == "0h_vs_2h",]
res2$dataset = factor(res2$dataset, levels = c("GSE90569", "GSE96538", "GSE94396", "GSE52260","GSE140244"))
p2 <- ggplot(res2, aes(x=dataset , y = gene , fill=log2FoldChange))+
theme_bw() +
geom_tile() +
scale_fill_distiller(palette = "RdBu", limits = c(-1.60, 1.60)) +
geom_text(data=res2, aes(x = dataset, y = gene, label = padj, angle = 0), size=2.0) +
scale_x_discrete(drop = FALSE) +
scale_y_discrete(drop = FALSE) +
theme(panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      text = element_text(size = 5),
      axis.text.x = element_text(angle = 0, size = 5),
      axis.text.y = element_text(angle = 0, size = 5),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 6, face="plain"),
      plot.tag = element_text(size = 8, face="plain"),
      legend.position="right",
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5),
      legend.title=element_text(size = 5)) +
labs(title = "0h vs 2h", tag = "B")

g<- ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="right")

# g <- arrangeGrob(p1, p2, ncol = 7, nrow = 1, layout_matrix= rbind(c(1,1,1,2,2,2,2)))
ggsave(paste0("time_conc_de_supplementary.png"), plot = g, units = "mm", width = 170, height = 60)



colnames(chaser_expr_d1) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')




t_cells <- read.csv('GSE140244_rnaseq_gene_counts.txt', sep = '\t', row.names = 1)
t_cells <- read.csv('GSE197067_HTSeq_counts.csv', row.names = 1)
