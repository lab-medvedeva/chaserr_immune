library(ggplot2)
library(GEOquery)
library(gridExtra)
library(stringr)
library(AnnotationHub)
library(biomaRt)
require("ensembldb")
library(ggpubr)

hub <- AnnotationHub()
query(hub, c("homo sapiens","ensdb"))
ensdb <- hub[["AH116291"]]
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl <- mart

get_tpm <- function(counts) 
{
  gene_lengths =  getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id', 'transcript_length','cds_length'), filters =  'ensembl_gene_id', values = rownames(counts), mart = ensembl, useCache = TRUE)
  gene_canonical_transcript =  getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','transcript_is_canonical'), filters = 'ensembl_gene_id', values = rownames(counts), mart = ensembl, useCache = TRUE)
  gene_canonical_transcript_subset = gene_canonical_transcript[!is.na(gene_canonical_transcript$transcript_is_canonical),]
  gene_lengths = merge(gene_canonical_transcript_subset, gene_lengths, by = c("ensembl_gene_id", "ensembl_transcript_id"))
  gene_lengths <- gene_lengths[!is.na(gene_lengths$transcript_length),]
  x <- counts[gene_lengths$ensembl_gene_id,]
  x <- x / gene_lengths$transcript_length
  counts_tpm <- t( t(x) * 1e6 / colSums(x) )
  return(counts_tpm)
}

# dataset GSE90569
treg <- read.csv('GSE90569_norm_counts_TPM_GRCh38.p13_NCBI.tsv', sep = '\t', row.names = 1)
treg_metadata <- read.csv('treg_metadata.csv', sep = '\t')
chaser_expr <- cbind(treg_metadata, t(treg['100507217',]), 'CHASERR')
chd2_expr <- cbind(treg_metadata, t(treg['1106',]), 'CHD2')
colnames(chaser_expr) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(chd2_expr) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
df_expr <- rbind(chaser_expr, chd2_expr)
Thp <- df_expr[df_expr$cell_type =='Thp',]
Thp_th0 <- df_expr
Thp_th0$cell_type <- 'Th0'
Thp_itreg <- df_expr
Thp_itreg$cell_type <- 'iTreg'
df_expr <- df_expr[df_expr$cell_type !='Thp',]
df_expr <- rbind(df_expr, Thp_itreg, Thp_th0)
df_expr$time <- paste0(df_expr$time, 'h')
df_expr$time <- factor(df_expr$time, levels = c('0h', '0.5h', '1h', '2h', '4h', '6h', '12h', '24h', '48h', '72h'))
p1 <- ggplot(df_expr, aes(x = time, y = expr, fill = gene)) + # , fill = celltype
  facet_grid(.~cell_type, scales = "fixed", space = "fixed") +
  theme_bw() +
  geom_boxplot() +
theme(
  strip.background = element_blank(),  
  strip.text.x = element_text(size = 6),
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  text = element_text(size = 6),
  axis.title.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.text.x = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  axis.title.y = element_text(size = 6),
  plot.title = element_text(size = 6, face="plain"),
  plot.tag = element_text(size = 8, face="plain"),
  legend.key.size = unit(0.3, "cm"),  # Key size
  legend.text = element_text(size = 6),  # Text labels
  legend.title = element_text(size = 6), # Title
  legend.position="right") +
  labs(title = "GSE90569", y = "Expression (TPM)", tag = "")

# dataset GSE96538
treg2 <- read.csv('GSE96538_norm_counts_TPM_GRCh38.p13_NCBI.tsv', sep = '\t', row.names = 1)
treg_metadata <- read.csv('treg2_metadata.csv', sep = '\t')
chaser_expr <- cbind(treg_metadata, t(treg2['100507217',]), 'CHASERR')
chd2_expr <- cbind(treg_metadata, t(treg2['1106',]), 'CHD2')
colnames(chaser_expr) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(chd2_expr) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
df_expr <- rbind(chaser_expr, chd2_expr)
Thp <- df_expr[df_expr$cell_type =='Thp',]
Thp_th0 <- df_expr
Thp_th0$cell_type <- 'Th0'
Thp_itreg <- df_expr
Thp_itreg$cell_type <- 'iTreg'
df_expr <- df_expr[df_expr$cell_type !='Thp',]
df_expr <- rbind(df_expr, Thp_itreg, Thp_th0)
df_expr$time <- paste0(df_expr$time, 'h')
df_expr$time <- factor(df_expr$time, levels = c('0h', '0.5h', '1h', '2h', '6h', '12h', '24h', '48h', '72h', '96h', '120h'))
p2 <- ggplot(df_expr, aes(x = time, y = expr, fill = gene)) + # , fill = celltype
  facet_grid(.~cell_type, scales = "fixed", space = "fixed") +
  theme_bw() +
  geom_boxplot() +
theme(
  strip.background = element_blank(),  
  strip.text.x = element_text(size = 6),
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  text = element_text(size = 6),
  axis.title.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.text.x = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  axis.title.y = element_text(size = 6),
  plot.title = element_text(size = 6, face="plain"),
  plot.tag = element_text(size = 8, face="plain"),
  legend.key.size = unit(0.3, "cm"),  # Key size
  legend.text = element_text(size = 6),  # Text labels
  legend.title = element_text(size = 6), # Title
  legend.position="right") +
  labs(title = "GSE96538", y = "Expression (TPM)", tag = "")

# dataset GSE94396
treg3 <- read.csv('GSE94396_norm_counts_TPM_GRCh38.p13_NCBI.tsv', sep = '\t', row.names = 1)
gse <- getGEO("GSE94396", GSEMatrix = TRUE)
gse_matrix <- pData(gse[[1]])
treg_metadata <- cbind(gse_matrix$geo_accession, gse_matrix$source_name_ch1, gse_matrix$'sample.time:ch1')
chaser_expr <- cbind(treg_metadata, t(treg3['100507217',]), 'CHASERR')
chd2_expr <- cbind(treg_metadata, t(treg3['1106',]), 'CHD2')
colnames(chaser_expr) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
colnames(chd2_expr) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
df_expr <- as.data.frame(rbind(chaser_expr, chd2_expr))
df_expr <- df_expr[df_expr$cell_type != "CD25high nTregs",]
df_expr$time <- paste0(df_expr$time, 'h')
df_expr$time <- factor(df_expr$time, levels = c("0h", "2h", "6h", "24h", "48h", "144h"))
df_expr$expr <- as.double(df_expr$expr)
p3 <- ggplot(df_expr, aes(x = time, y = expr, fill = gene)) + # , fill = celltype
  theme_bw() +
  geom_boxplot() +
theme(
  strip.background = element_blank(),  
  strip.text.x = element_text(size = 6),
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  text = element_text(size = 6),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.x = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  axis.title.y = element_text(size = 6),
  plot.title = element_text(size = 6, face="plain"),
  plot.tag = element_text(size = 8, face="plain"),
  legend.key.size = unit(0.3, "cm"),  # Key size
  legend.text = element_text(size = 6),  # Text labels
  legend.title = element_text(size = 6), # Title
  legend.position="right")  +
  labs(title = "GSE94396", y = "Expression (TPM)", tag = "")

# dataset GSE52260
th17 <- read.csv('GSE52260_norm_counts_TPM_GRCh38.p13_NCBI.tsv', sep = '\t', row.names = 1)
th17_metadata <- read.csv('th17_metadata_human.csv', sep = '\t')
chaser_expr <- cbind(th17_metadata, t(th17['100507217',]), 'CHASERR')
chd2_expr <- cbind(th17_metadata, t(th17['1106',]), 'CHD2')
colnames(chaser_expr) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(chd2_expr) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
df_expr <- rbind(chaser_expr, chd2_expr)
Thp <- df_expr[df_expr$cell_type =='Thp',]
Thp_th0 <- df_expr
Thp_th0$cell_type <- 'Th0'
Thp_th17 <- df_expr
Thp_th17$cell_type <- 'Th17'
df_expr <- df_expr[df_expr$cell_type !='Thp',]
df_expr <- rbind(df_expr, Thp_th17, Thp_th0)
df_expr$time <- factor(df_expr$time, levels = c('0h', '0.5h', '1h', '2h', '4h', '6h', '12h', '24h', '48h', '72h'))
p4 <- ggplot(df_expr, aes(x = time, y = expr, fill = gene)) + # , fill = celltype
  facet_grid(.~cell_type, scales = "fixed", space = "fixed") +
  theme_bw() +
  geom_boxplot() +
theme(
  strip.background = element_blank(),  
  strip.text.x = element_text(size = 6),
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  text = element_text(size = 6),
  axis.title.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.text.x = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  axis.title.y = element_text(size = 6),
  plot.title = element_text(size = 6, face="plain"),
  plot.tag = element_text(size = 8, face="plain"),
  legend.key.size = unit(0.3, "cm"),  # Key size
  legend.text = element_text(size = 6),  # Text labels
  legend.title = element_text(size = 6), # Title
  legend.position="right") +
  labs(title = "GSE52260", y = "Expression (TPM)", tag = "")

# dataset GSE140244
t_cells <- read.csv('GSE140244_rnaseq_gene_counts.txt', sep = '\t', row.names = 1)
counts_tpm <- get_tpm(t_cells)
counts_chaserr_tpm <- counts_tpm['ENSG00000272888',]
counts_chd2_tpm <- counts_tpm['ENSG00000173575',]
metadata <- read.csv('GSE140244_rnaseq_meta_data.txt', sep = '\t')
chd2_expr <- cbind(metadata$ExpressionMatrix_SampleID, metadata$Time_point, as.vector(counts_chd2_tpm), 'CHD2') # CHD2
colnames(chd2_expr) <- c('sample', 'time', 'expr', 'gene')
chaser_expr <- cbind(metadata$ExpressionMatrix_SampleID, metadata$Time_point, as.vector(counts_chaserr_tpm), 'CHASERR') # CHASERR
colnames(chaser_expr) <- c('sample', 'time', 'expr', 'gene')
df_expr <- as.data.frame(rbind(chaser_expr, chd2_expr))
df_expr$time <- paste0(df_expr$time, 'h')
df_expr$time <- factor(df_expr$time, levels = c("0h", "2h", "4h", "8h", "12h", "24h", "48h", "72h"))
df_expr$expr <-as.double(df_expr$expr)
p5 <- ggplot(df_expr, aes(x = time, y = expr, fill = gene)) + # , fill = celltype
  theme_bw() +
  geom_boxplot()  +
theme(
  strip.background = element_blank(),  
  strip.text.x = element_text(size = 6),
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  text = element_text(size = 6),
  axis.title.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.text.x = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  axis.title.y = element_text(size = 6),
  plot.title = element_text(size = 6, face="plain"),
  plot.tag = element_text(size = 8, face="plain"),
  legend.key.size = unit(0.3, "cm"),  # Key size
  legend.text = element_text(size = 6),  # Text labels
  legend.title = element_text(size = 6), # Title
  legend.position="right") +
  labs(title = "GSE140244", y = "Expression (TPM)", tag = "")

# dataset GSE197067
t_cells <- read.csv('GSE197067_HTSeq_counts.csv', row.names = 1)
rownames(t_cells) <- str_split_fixed(rownames(t_cells), "\\.", 2)[,1]
counts_tpm <- get_tpm(t_cells)
pan_t_metadata <- read.csv('pan_t_metadata.csv', sep = '\t')
pan_t_metadata <- pan_t_metadata[pan_t_metadata$cell_type == 'act',]
counts_tpm <- counts_tpm[pan_t_metadata$cell_type == 'act',]
counts_chaserr_tpm <- counts_tpm['ENSG00000272888',]
counts_chd2_tpm <- counts_tpm['ENSG00000173575',]
chaser_expr <- cbind(pan_t_metadata$sample, pan_t_metadata$cell_type, pan_t_metadata$time, as.vector(counts_chaserr_tpm), 'CHASERR') # CHASERR
colnames(chaser_expr) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
chd2_expr <- cbind(pan_t_metadata$sample, pan_t_metadata$cell_type, pan_t_metadata$time, as.vector(counts_chd2_tpm), 'CHD2') # CHD2
colnames(chd2_expr) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
df_expr <- as.data.frame(rbind(chaser_expr, chd2_expr))
df_expr$time <- factor(df_expr$time, levels = c("0h", "6h", "12h", "24h", "48h", "72h"))
df_expr$expr <-as.double(df_expr$expr)
p6 <- ggplot(df_expr, aes(x = time, y = expr, fill = gene)) + # , fill = celltype
  theme_bw() +
  geom_boxplot() +
theme(
  strip.background = element_blank(),  
  strip.text.x = element_text(size = 6),
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  text = element_text(size = 6),
  axis.title.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.text.x = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  axis.title.y = element_text(size = 6),
  plot.title = element_text(size = 6, face="plain"),
  plot.tag = element_text(size = 8, face="plain"),
  legend.key.size = unit(0.3, "cm"),  # Key size
  legend.text = element_text(size = 6),  # Text labels
  legend.title = element_text(size = 6), # Title
  legend.position="right") +
  labs(title = "GSE197067", y = "Expression (TPM)", tag = "")


h <- ggarrange(p1, p2, p3, p4, p5, p6, align='v',
          nrow=6, ncol = 1, common.legend = TRUE, legend="bottom")


# g <- arrangeGrob(p1, p2, p3, p4, p5, p6, ncol = 1, nrow = 6, layout_matrix= rbind(c(1), c(2), c(3), c(4), c(5), c(6)))
ggsave(paste0("time_activation_supplementary.png"), plot = h, units = "mm", width = 170, height = 225)

 
