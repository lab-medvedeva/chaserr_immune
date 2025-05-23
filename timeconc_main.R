library(ggplot2)
library(GEOquery)
library(gridExtra)
library(stringr)
library(AnnotationHub)
library(biomaRt)
require("ensembldb")
library(ggpubr)

# B. T cell activation expression

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
chaser_expr_d1 <- cbind(treg_metadata, t(treg['100507217',]), 'CHASERR')
chd2_expr_d1 <- cbind(treg_metadata, t(treg['1106',]), 'CHD2')
colnames(chaser_expr_d1) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(chd2_expr_d1) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
chaser_expr_d1$rep <- NULL
chd2_expr_d1$rep <- NULL
chaser_expr_d1$expr <- scale(as.double(chaser_expr_d1$expr))
chd2_expr_d1$expr <- scale(as.double(chd2_expr_d1$expr))

# dataset GSE96538
treg2 <- read.csv('GSE96538_norm_counts_TPM_GRCh38.p13_NCBI.tsv', sep = '\t', row.names = 1)
treg_metadata <- read.csv('treg2_metadata.csv', sep = '\t')
chaser_expr_d2 <- cbind(treg_metadata, t(treg2['100507217',]), 'CHASERR')
chd2_expr_d2 <- cbind(treg_metadata, t(treg2['1106',]), 'CHD2')
colnames(chaser_expr_d2) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(chd2_expr_d2) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
chaser_expr_d2$rep <- NULL
chd2_expr_d2$rep <- NULL
chaser_expr_d2$expr <- scale(as.double(chaser_expr_d2$expr))
chd2_expr_d2$expr <- scale(as.double(chd2_expr_d2$expr))

# dataset GSE94396
treg3 <- read.csv('GSE94396_norm_counts_TPM_GRCh38.p13_NCBI.tsv', sep = '\t', row.names = 1)
gse <- getGEO("GSE94396", GSEMatrix = TRUE)
gse_matrix <- pData(gse[[1]])
treg_metadata <- read.csv('GSE94396_series_matrix.txt', sep = '\t')
chaser_expr_d3 <- cbind(treg_metadata, t(treg3['100507217',]), 'CHASERR')
chd2_expr_d3 <- cbind(treg_metadata, t(treg3['1106',]), 'CHD2')
colnames(chaser_expr_d3) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
colnames(chd2_expr_d3) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
chaser_expr_d3 <- as.data.frame(chaser_expr_d3)
chd2_expr_d3 <- as.data.frame(chd2_expr_d3)
chaser_expr_d3[chaser_expr_d3$cell_type == "CD4+ T cells", 'cell_type'] = 'Th0'
chaser_expr_d3[chaser_expr_d3$cell_type == "CD25high nTregs", 'cell_type'] = 'Treg'
chd2_expr_d3[chd2_expr_d3$cell_type == "CD4+ T cells", 'cell_type'] = 'Th0'
chd2_expr_d3[chd2_expr_d3$cell_type == "CD25high nTregs", 'cell_type'] = 'Treg'
chaser_expr_d3$expr <- scale(as.double(chaser_expr_d3$expr))
chd2_expr_d3$expr <- scale(as.double(chd2_expr_d3$expr))

# dataset GSE52260
th17 <- read.csv('GSE52260_norm_counts_TPM_GRCh38.p13_NCBI.tsv', sep = '\t', row.names = 1)
th17_metadata <- read.csv('th17_metadata_human.csv', sep = '\t')
chaser_expr_d4 <- cbind(th17_metadata, t(th17['100507217',]), 'CHASERR')
chd2_expr_d4 <- cbind(th17_metadata, t(th17['1106',]), 'CHD2')
colnames(chaser_expr_d4) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(chd2_expr_d4) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
chaser_expr_d4$rep <- NULL
chd2_expr_d4$rep <-NULL

chaser_expr_d4[chaser_expr_d4$time == '0h', 'time'] = '0'
chaser_expr_d4[chaser_expr_d4$time == '0.5h', 'time'] = '0.5'
chaser_expr_d4[chaser_expr_d4$time == '1h', 'time'] = '1'
chaser_expr_d4[chaser_expr_d4$time == '2h', 'time'] = '2'
chaser_expr_d4[chaser_expr_d4$time == '4h', 'time'] = '4'
chaser_expr_d4[chaser_expr_d4$time == '6h', 'time'] = '6'
chaser_expr_d4[chaser_expr_d4$time == '12h', 'time'] = '12'
chaser_expr_d4[chaser_expr_d4$time == '24h', 'time'] = '24'
chaser_expr_d4[chaser_expr_d4$time == '48h', 'time'] = '48'
chaser_expr_d4[chaser_expr_d4$time == '72h', 'time'] = '72'

chd2_expr_d4[chd2_expr_d4$time == '0h', 'time'] = '0'
chd2_expr_d4[chd2_expr_d4$time == '0.5h', 'time'] = '0.5'
chd2_expr_d4[chd2_expr_d4$time == '1h', 'time'] = '1'
chd2_expr_d4[chd2_expr_d4$time == '2h', 'time'] = '2'
chd2_expr_d4[chd2_expr_d4$time == '4h', 'time'] = '4'
chd2_expr_d4[chd2_expr_d4$time == '6h', 'time'] = '6'
chd2_expr_d4[chd2_expr_d4$time == '12h', 'time'] = '12'
chd2_expr_d4[chd2_expr_d4$time == '24h', 'time'] = '24'
chd2_expr_d4[chd2_expr_d4$time == '48h', 'time'] = '48'
chd2_expr_d4[chd2_expr_d4$time == '72h', 'time'] = '72'

chaser_expr_d4$expr <- scale(as.double(chaser_expr_d4$expr))
chd2_expr_d4$expr <- scale(as.double(chd2_expr_d4$expr))

# dataset GSE140244
t_cells <- read.csv('GSE140244_rnaseq_gene_counts.txt', sep = '\t', row.names = 1)
counts_tpm <- get_tpm(t_cells)
counts_chaserr_tpm <- counts_tpm['ENSG00000272888',]
counts_chd2_tpm <- counts_tpm['ENSG00000173575',]
metadata <- read.csv('GSE140244_rnaseq_meta_data.txt', sep = '\t')
chd2_expr_d5 <- cbind(metadata$ExpressionMatrix_SampleID, 'pan T', metadata$Time_point, as.vector(counts_chd2_tpm), 'CHD2') # CHD2
colnames(chd2_expr_d5) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
chaser_expr_d5 <- cbind(metadata$ExpressionMatrix_SampleID, 'pan T', metadata$Time_point, as.vector(counts_chaserr_tpm), 'CHASERR') # CHASERR
colnames(chaser_expr_d5) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
chaser_expr_d5 <- as.data.frame(chaser_expr_d5)
chd2_expr_d5 <- as.data.frame(chd2_expr_d5)
chaser_expr_d5$expr <- scale(as.double(chaser_expr_d5$expr))
chd2_expr_d5$expr <- scale(as.double(chd2_expr_d5$expr))

# dataset GSE197067
t_cells <- read.csv('GSE197067_HTSeq_counts.csv', row.names = 1)
rownames(t_cells) <- str_split_fixed(rownames(t_cells), "\\.", 2)[,1]
counts_tpm <- get_tpm(t_cells)
pan_t_metadata <- read.csv('pan_t_metadata.csv', sep = '\t')
pan_t_metadata <- pan_t_metadata[pan_t_metadata$cell_type == 'act',]
counts_tpm <- counts_tpm[pan_t_metadata$cell_type == 'act',]
counts_chaserr_tpm <- counts_tpm['ENSG00000272888',]
counts_chd2_tpm <- counts_tpm['ENSG00000173575',]

chaser_expr_d6 <- cbind(pan_t_metadata$sample, pan_t_metadata$cell_type, pan_t_metadata$time, as.vector(counts_chaserr_tpm), 'CHASERR') # CHASERR
colnames(chaser_expr_d6) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
chd2_expr_d6 <- cbind(pan_t_metadata$sample, pan_t_metadata$cell_type, pan_t_metadata$time, as.vector(counts_chd2_tpm), 'CHD2') # CHD2
colnames(chd2_expr_d6) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
chaser_expr_d6 <- as.data.frame(chaser_expr_d6)
chd2_expr_d6 <- as.data.frame(chd2_expr_d6)

chaser_expr_d6[chaser_expr_d6$time == '0h', 'time'] = '0'
chaser_expr_d6[chaser_expr_d6$time == '6h', 'time'] = '6'
chaser_expr_d6[chaser_expr_d6$time == '12h', 'time'] = '12'
chaser_expr_d6[chaser_expr_d6$time == '24h', 'time'] = '24'
chaser_expr_d6[chaser_expr_d6$time == '48h', 'time'] = '48'
chaser_expr_d6[chaser_expr_d6$time == '72h', 'time'] = '72'

chd2_expr_d6[chd2_expr_d6$time == '0h', 'time'] = '0'
chd2_expr_d6[chd2_expr_d6$time == '6h', 'time'] = '6'
chd2_expr_d6[chd2_expr_d6$time == '12h', 'time'] = '12'
chd2_expr_d6[chd2_expr_d6$time == '24h', 'time'] = '24'
chd2_expr_d6[chd2_expr_d6$time == '48h', 'time'] = '48'
chd2_expr_d6[chd2_expr_d6$time == '72h', 'time'] = '72'

chaser_expr_d6$expr <- scale(as.double(chaser_expr_d6$expr))
chd2_expr_d6$expr <- scale(as.double(chd2_expr_d6$expr))

# concat data
chaser_expr <- rbind(chaser_expr_d1, chaser_expr_d2, chaser_expr_d3, chaser_expr_d4, chaser_expr_d5, chaser_expr_d6)#)
chaser_expr$expr <- as.double(chaser_expr$expr)
chaser_expr$time <- factor(as.character(chaser_expr$time), levels = c( "0", "0.5", "1", "2", "4", "6", "8", "12", "24", "48", "72", "96", "120", "144"))
chd2_expr <- rbind(chd2_expr_d1, chd2_expr_d2, chd2_expr_d3, chd2_expr_d4, chd2_expr_d5, chaser_expr_d6)
chd2_expr$expr <- as.double(chd2_expr$expr)
chd2_expr$time <- factor(as.character(chd2_expr$time), levels = c( "0", "0.5", "1", "2", "4", "6", "8", "12", "24", "48", "72", "96", "120", "144"))
all_expr <- rbind(chaser_expr, chd2_expr)

p6 <- ggplot(all_expr, aes(x = time, y = expr, fill = gene)) + # , fill = celltype
  facet_grid(.~gene, scales = "fixed", space = "fixed") +
  theme_bw() +
  geom_boxplot() +
  theme(panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  text = element_text(size = 6),
  axis.text.x = element_text(angle = 90, size = 6),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  legend.title = element_blank(), # Title
  plot.title = element_text(size = 6, face = "plain"),
  plot.tag = element_text(size = 8, face = "plain")) +
  geom_smooth(method="loess", aes(group = 1)) +
  labs(x = "Time (hours)" , y = "Expression (TPM)", title = "", tag = "B")

# B. T cell activation and FANTOM6 knockdown signature
# T cell activation differentially expressed genes
# https://doi.org/10.1186/s13059-023-03120-7
# https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-023-03120-7/MediaObjects/13059_2023_3120_MOESM2_ESM.html

modules_activation <- read.csv('T_cell_activation.csv', sep = ',')
rownames(modules_activation) <- modules_activation$Ensembl

# FANTOM6 signature
# https://fantom.gsc.riken.jp/6/datafiles/Core_FANTOM6/RELEASE_latest/analysis/DEGs/01_combined/DESeq2_genes_ASO_all.tsv.bz2

aso_df <- read.csv("../DESeq2_genes_ASO_all.tsv", sep = '\t')
aso_df_fdr <- dplyr::filter(aso_df, fdr < 0.05)
aso_df_07 <- dplyr::filter(aso_df_fdr, perturb_id == 'ASO_G0272888_AD_07')
aso_df_10 <- dplyr::filter(aso_df_fdr, perturb_id == 'ASO_G0272888_AD_10')
inter_07 <- intersect(aso_df_07$geneID, modules_activation$Ensembl)
rownames(aso_df_07) <- aso_df_07$geneID
aso_df_07_module <- aso_df_07[inter_07,]
aso_df_07_module$aso <- 'ASO_07'
aso_df_07_module$Effect <- modules_activation[inter_07, 'Effect']
inter_10 <- intersect(aso_df_10$geneID, modules_activation$Ensembl)
rownames(aso_df_10) <- aso_df_10$geneID
aso_df_10_module <- aso_df_10[inter_10,]
aso_df_10_module$aso <- 'ASO_10'
aso_df_10_module$Effect <- modules_activation[inter_10, 'Effect']
modules_activation_inter <- rbind(aso_df_07_module, aso_df_10_module)
modules_activation_inter$to_include <- sign(modules_activation_inter$log2FC) == sign(modules_activation_inter$Effect)
write.table(modules_activation_inter, file = paste0("modules_activation_inter.csv"), sep="\t", col.names=TRUE, quote = FALSE)
modules_activation_inter <- modules_activation_inter[order(modules_activation_inter$to_include), ]
modules_activation_inter <- modules_activation_inter[modules_activation_inter$geneSymbol != "AC013394.2", ]
modules_activation_inter1 <- modules_activation_inter
modules_activation_inter1$type <- 'FANTOM6 logfc'

p1 <- ggplot(modules_activation_inter1, aes(x = geneSymbol, y = type, fill = log2FC)) +
facet_grid(. ~ aso, scales = "free",  space = "free") +
   geom_tile(width = 1) +
   theme_bw() +
scale_fill_distiller(palette = "RdBu", limits = c(-1, 1)) +
scale_x_discrete(drop = TRUE) +
scale_y_discrete(drop = FALSE) +
theme(
  strip.background = element_blank(),  
  strip.text.x = element_text(size = 6),
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  text = element_text(size = 6),
  axis.text.x = element_blank(),
  axis.title.x = element_blank(),
  axis.ticks.x=element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank(),
  plot.title = element_text(size = 6, face = "plain"),
  plot.tag = element_text(size = 8, face = "plain"),
  legend.key.size = unit(0.3, "cm"),
  legend.text = element_text(size = 6),
  legend.title = element_text(size = 6),
  legend.position="right") +
  labs(title = "", tag = "C")

legend1 <- cowplot::get_legend(p1)

modules_activation_inter2 <- modules_activation_inter
modules_activation_inter2$type <- 'Activation Effect'

p2 <- ggplot(modules_activation_inter2, aes(x = geneSymbol, y = type, fill = Effect )) +
facet_grid(. ~ aso, scales = "free",  space = "free") +
   geom_tile(width = 1) +
   theme_bw() +
scale_fill_distiller(palette = "RdBu", limits = c(-8, 8)) +
scale_x_discrete(drop = TRUE) +
scale_y_discrete(drop = FALSE) +
theme(panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  text = element_text(size = 6),
  axis.text.x = element_text(angle = 0, size = 5),
  axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
        legend.key.size = unit(0.3, "cm"),  # Key size
  legend.text = element_text(size = 6),  # Text labels
  legend.title = element_text(size = 6), # Title
      legend.position="right") 

legend2 <- cowplot::get_legend(p2)

p <- ggarrange(p1, p2, align='v',
          nrow=2, ncol = 1, common.legend = TRUE, legend="none")

g <- arrangeGrob(p6, p, legend1, legend2,  ncol = 10, nrow = 4 , 
  layout_matrix= rbind(c(1,1,1,1,1,1,1,1,1,1), c(1,1,1,1,1,1,1,1,1,1), c(2,2,2,2,2,2,2,2,2,3), c(2,2,2,2,2,2,2,2,2,4)))

saveRDS(g, "time_activation_main.rds")
