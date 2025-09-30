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
rictor_expr_d1 <- cbind(treg_metadata, t(treg['253260',]), 'RICTOR')
rbl2_expr_d1 <- cbind(treg_metadata, t(treg['5934',]), 'RBL2')
usp30_expr_d1 <- cbind(treg_metadata, t(treg['84749',]), 'USP30')
hn1_expr_d1 <- cbind(treg_metadata, t(treg['51155',]), 'HN1')
osbpl3_expr_d1 <- cbind(treg_metadata, t(treg['26031',]), 'OSBPL3')

colnames(rictor_expr_d1) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(rbl2_expr_d1) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(usp30_expr_d1) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(hn1_expr_d1) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(osbpl3_expr_d1) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')

rictor_expr_d1$rep <- NULL
rbl2_expr_d1$rep <- NULL
usp30_expr_d1$rep <- NULL
hn1_expr_d1$rep <- NULL
osbpl3_expr_d1$rep <- NULL

rictor_expr_d1$expr <- scale(as.double(rictor_expr_d1$expr))
rbl2_expr_d1$expr <- scale(as.double(rbl2_expr_d1$expr))
usp30_expr_d1$expr <- scale(as.double(usp30_expr_d1$expr))
hn1_expr_d1$expr <- scale(as.double(hn1_expr_d1$expr))
osbpl3_expr_d1$expr <- scale(as.double(osbpl3_expr_d1$expr))

# dataset GSE96538
treg2 <- read.csv('GSE96538_norm_counts_TPM_GRCh38.p13_NCBI.tsv', sep = '\t', row.names = 1)
treg_metadata <- read.csv('treg2_metadata.csv', sep = '\t')
rictor_expr_d2 <- cbind(treg_metadata, t(treg2['253260',]), 'RICTOR')
rbl2_expr_d2 <- cbind(treg_metadata, t(treg2['5934',]), 'RBL2')
usp30_expr_d2 <- cbind(treg_metadata, t(treg2['84749',]), 'USP30')
hn1_expr_d2 <- cbind(treg_metadata, t(treg2['51155',]), 'HN1')
osbpl3_expr_d2 <- cbind(treg_metadata, t(treg2['26031',]), 'OSBPL3')

colnames(rictor_expr_d2) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(rbl2_expr_d2) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(usp30_expr_d2) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(hn1_expr_d2) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(osbpl3_expr_d2) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')

rictor_expr_d2$rep <- NULL
rbl2_expr_d2$rep <- NULL
usp30_expr_d2$rep <- NULL
hn1_expr_d2$rep <- NULL
osbpl3_expr_d2$rep <- NULL

rictor_expr_d2$expr <- scale(as.double(rictor_expr_d2$expr))
rbl2_expr_d2$expr <- scale(as.double(rbl2_expr_d2$expr))
usp30_expr_d2$expr <- scale(as.double(usp30_expr_d2$expr))
hn1_expr_d2$expr <- scale(as.double(hn1_expr_d2$expr))
osbpl3_expr_d2$expr <- scale(as.double(osbpl3_expr_d2$expr))

# dataset GSE94396
treg3 <- read.csv('GSE94396_norm_counts_TPM_GRCh38.p13_NCBI.tsv', sep = '\t', row.names = 1)
gse <- getGEO("GSE94396", GSEMatrix = TRUE)
gse_matrix <- pData(gse[[1]])
treg_metadata <- read.csv('GSE94396_series_matrix.txt', sep = '\t')

rictor_expr_d3 <- cbind(treg_metadata, t(treg3['253260',]), 'RICTOR')
rbl2_expr_d3 <- cbind(treg_metadata, t(treg3['5934',]), 'RBL2')
usp30_expr_d3 <- cbind(treg_metadata, t(treg3['84749',]), 'USP30')
hn1_expr_d3 <- cbind(treg_metadata, t(treg3['51155',]), 'HN1')
osbpl3_expr_d3 <- cbind(treg_metadata, t(treg3['26031',]), 'OSBPL3')

colnames(rictor_expr_d3) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
colnames(rbl2_expr_d3) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
colnames(usp30_expr_d3) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
colnames(hn1_expr_d3) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
colnames(osbpl3_expr_d3) <- c('sample', 'cell_type', 'time', 'expr', 'gene')

rictor_expr_d3 <- as.data.frame(rictor_expr_d3)
rbl2_expr_d3 <- as.data.frame(rbl2_expr_d3)
usp30_expr_d3 <- as.data.frame(usp30_expr_d3)
hn1_expr_d3 <- as.data.frame(hn1_expr_d3)
osbpl3_expr_d3 <- as.data.frame(osbpl3_expr_d3)

rictor_expr_d3[rictor_expr_d3$cell_type == "CD4+ T cells", 'cell_type'] = 'Th0'
rbl2_expr_d3[rbl2_expr_d3$cell_type == "CD4+ T cells", 'cell_type'] = 'Th0'
usp30_expr_d3[usp30_expr_d3$cell_type == "CD4+ T cells", 'cell_type'] = 'Th0'
hn1_expr_d3[hn1_expr_d3$cell_type == "CD4+ T cells", 'cell_type'] = 'Th0'
osbpl3_expr_d3[osbpl3_expr_d3$cell_type == "CD4+ T cells", 'cell_type'] = 'Th0'

rictor_expr_d3[rictor_expr_d3$cell_type == "CD25high nTregs", 'cell_type'] = 'Treg'
rbl2_expr_d3[rbl2_expr_d3$cell_type == "CD25high nTregs", 'cell_type'] = 'Treg'
usp30_expr_d3[usp30_expr_d3$cell_type == "CD25high nTregs", 'cell_type'] = 'Treg'
hn1_expr_d3[hn1_expr_d3$cell_type == "CD25high nTregs", 'cell_type'] = 'Treg'
osbpl3_expr_d3[osbpl3_expr_d3$cell_type == "CD25high nTregs", 'cell_type'] = 'Treg'

rictor_expr_d3$expr <- scale(as.double(rictor_expr_d3$expr))
rbl2_expr_d3$expr <- scale(as.double(rbl2_expr_d3$expr))
usp30_expr_d3$expr <- scale(as.double(usp30_expr_d3$expr))
hn1_expr_d3$expr <- scale(as.double(hn1_expr_d3$expr))
osbpl3_expr_d3$expr <- scale(as.double(osbpl3_expr_d3$expr))

# dataset GSE52260
th17 <- read.csv('GSE52260_norm_counts_TPM_GRCh38.p13_NCBI.tsv', sep = '\t', row.names = 1)
th17_metadata <- read.csv('th17_metadata_human.csv', sep = '\t')

rictor_expr_d4 <- cbind(th17_metadata, t(th17['253260',]), 'RICTOR')
rbl2_expr_d4 <- cbind(th17_metadata, t(th17['5934',]), 'RBL2')
usp30_expr_d4 <- cbind(th17_metadata, t(th17['84749',]), 'USP30')
hn1_expr_d4 <- cbind(th17_metadata, t(th17['51155',]), 'HN1')
osbpl3_expr_d4 <- cbind(th17_metadata, t(th17['26031',]), 'OSBPL3')

colnames(rictor_expr_d4) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(rbl2_expr_d4) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(usp30_expr_d4) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(hn1_expr_d4) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')
colnames(osbpl3_expr_d4) <- c('sample', 'cell_type', 'time', 'rep', 'expr', 'gene')

rictor_expr_d4$rep <- NULL
rbl2_expr_d4$rep <-NULL
usp30_expr_d4$rep <- NULL
hn1_expr_d4$rep <-NULL
osbpl3_expr_d4$rep <- NULL

rictor_expr_d4[rictor_expr_d4$time == '0h', 'time'] = '0'
rictor_expr_d4[rictor_expr_d4$time == '0.5h', 'time'] = '0.5'
rictor_expr_d4[rictor_expr_d4$time == '1h', 'time'] = '1'
rictor_expr_d4[rictor_expr_d4$time == '2h', 'time'] = '2'
rictor_expr_d4[rictor_expr_d4$time == '4h', 'time'] = '4'
rictor_expr_d4[rictor_expr_d4$time == '6h', 'time'] = '6'
rictor_expr_d4[rictor_expr_d4$time == '12h', 'time'] = '12'
rictor_expr_d4[rictor_expr_d4$time == '24h', 'time'] = '24'
rictor_expr_d4[rictor_expr_d4$time == '48h', 'time'] = '48'
rictor_expr_d4[rictor_expr_d4$time == '72h', 'time'] = '72'

rbl2_expr_d4[rbl2_expr_d4$time == '0h', 'time'] = '0'
rbl2_expr_d4[rbl2_expr_d4$time == '0.5h', 'time'] = '0.5'
rbl2_expr_d4[rbl2_expr_d4$time == '1h', 'time'] = '1'
rbl2_expr_d4[rbl2_expr_d4$time == '2h', 'time'] = '2'
rbl2_expr_d4[rbl2_expr_d4$time == '4h', 'time'] = '4'
rbl2_expr_d4[rbl2_expr_d4$time == '6h', 'time'] = '6'
rbl2_expr_d4[rbl2_expr_d4$time == '12h', 'time'] = '12'
rbl2_expr_d4[rbl2_expr_d4$time == '24h', 'time'] = '24'
rbl2_expr_d4[rbl2_expr_d4$time == '48h', 'time'] = '48'
rbl2_expr_d4[rbl2_expr_d4$time == '72h', 'time'] = '72'

usp30_expr_d4[usp30_expr_d4$time == '0h', 'time'] = '0'
usp30_expr_d4[usp30_expr_d4$time == '0.5h', 'time'] = '0.5'
usp30_expr_d4[usp30_expr_d4$time == '1h', 'time'] = '1'
usp30_expr_d4[usp30_expr_d4$time == '2h', 'time'] = '2'
usp30_expr_d4[usp30_expr_d4$time == '4h', 'time'] = '4'
usp30_expr_d4[usp30_expr_d4$time == '6h', 'time'] = '6'
usp30_expr_d4[usp30_expr_d4$time == '12h', 'time'] = '12'
usp30_expr_d4[usp30_expr_d4$time == '24h', 'time'] = '24'
usp30_expr_d4[usp30_expr_d4$time == '48h', 'time'] = '48'
usp30_expr_d4[usp30_expr_d4$time == '72h', 'time'] = '72'

hn1_expr_d4[hn1_expr_d4$time == '0h', 'time'] = '0'
hn1_expr_d4[hn1_expr_d4$time == '0.5h', 'time'] = '0.5'
hn1_expr_d4[hn1_expr_d4$time == '1h', 'time'] = '1'
hn1_expr_d4[hn1_expr_d4$time == '2h', 'time'] = '2'
hn1_expr_d4[hn1_expr_d4$time == '4h', 'time'] = '4'
hn1_expr_d4[hn1_expr_d4$time == '6h', 'time'] = '6'
hn1_expr_d4[hn1_expr_d4$time == '12h', 'time'] = '12'
hn1_expr_d4[hn1_expr_d4$time == '24h', 'time'] = '24'
hn1_expr_d4[hn1_expr_d4$time == '48h', 'time'] = '48'
hn1_expr_d4[hn1_expr_d4$time == '72h', 'time'] = '72'

osbpl3_expr_d4[osbpl3_expr_d4$time == '0h', 'time'] = '0'
osbpl3_expr_d4[osbpl3_expr_d4$time == '0.5h', 'time'] = '0.5'
osbpl3_expr_d4[osbpl3_expr_d4$time == '1h', 'time'] = '1'
osbpl3_expr_d4[osbpl3_expr_d4$time == '2h', 'time'] = '2'
osbpl3_expr_d4[osbpl3_expr_d4$time == '4h', 'time'] = '4'
osbpl3_expr_d4[osbpl3_expr_d4$time == '6h', 'time'] = '6'
osbpl3_expr_d4[osbpl3_expr_d4$time == '12h', 'time'] = '12'
osbpl3_expr_d4[osbpl3_expr_d4$time == '24h', 'time'] = '24'
osbpl3_expr_d4[osbpl3_expr_d4$time == '48h', 'time'] = '48'
osbpl3_expr_d4[osbpl3_expr_d4$time == '72h', 'time'] = '72'

rictor_expr_d4$expr <- scale(as.double(rictor_expr_d4$expr))
rbl2_expr_d4$expr <- scale(as.double(rbl2_expr_d4$expr))
usp30_expr_d4$expr <- scale(as.double(usp30_expr_d4$expr))
hn1_expr_d4$expr <- scale(as.double(hn1_expr_d4$expr))
osbpl3_expr_d4$expr <- scale(as.double(osbpl3_expr_d4$expr))

# dataset GSE140244
t_cells <- read.csv('GSE140244_rnaseq_gene_counts.txt', sep = '\t', row.names = 1)
counts_tpm <- get_tpm(t_cells)

counts_rictor_tpm <- counts_tpm['ENSG00000164327',]
counts_rbl2_tpm <- counts_tpm['ENSG00000103479',]
counts_usp30_tpm <- counts_tpm['ENSG00000135093',]
counts_hn1_tpm <- counts_tpm['ENSG00000189159',]
counts_osbpl3_tpm <- counts_tpm['ENSG00000070882',]

metadata <- read.csv('GSE140244_rnaseq_meta_data.txt', sep = '\t')

rictor_expr_d5 <- cbind(metadata$ExpressionMatrix_SampleID, 'pan T', metadata$Time_point, as.vector(counts_rictor_tpm), 'RICTOR')
rbl2_expr_d5 <- cbind(metadata$ExpressionMatrix_SampleID, 'pan T', metadata$Time_point, as.vector(counts_rbl2_tpm), 'RBL2')
usp30_expr_d5 <- cbind(metadata$ExpressionMatrix_SampleID, 'pan T', metadata$Time_point, as.vector(counts_usp30_tpm), 'USP30')
hn1_expr_d5 <- cbind(metadata$ExpressionMatrix_SampleID, 'pan T', metadata$Time_point, as.vector(counts_hn1_tpm), 'HN1')
osbpl3_expr_d5 <- cbind(metadata$ExpressionMatrix_SampleID, 'pan T', metadata$Time_point, as.vector(counts_osbpl3_tpm), 'OSBPL3')

colnames(rictor_expr_d5) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
colnames(rbl2_expr_d5) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
colnames(usp30_expr_d5) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
colnames(hn1_expr_d5) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
colnames(osbpl3_expr_d5) <- c('sample', 'cell_type', 'time', 'expr', 'gene')

rictor_expr_d5 <- as.data.frame(rictor_expr_d5)
rbl2_expr_d5 <- as.data.frame(rbl2_expr_d5)
usp30_expr_d5 <- as.data.frame(usp30_expr_d5)
hn1_expr_d5 <- as.data.frame(hn1_expr_d5)
osbpl3_expr_d5 <- as.data.frame(osbpl3_expr_d5)

rictor_expr_d5$expr <- scale(as.double(rictor_expr_d5$expr))
rbl2_expr_d5$expr <- scale(as.double(rbl2_expr_d5$expr))
usp30_expr_d5$expr <- scale(as.double(usp30_expr_d5$expr))
hn1_expr_d5$expr <- scale(as.double(hn1_expr_d5$expr))
osbpl3_expr_d5$expr <- scale(as.double(osbpl3_expr_d5$expr))

# dataset GSE197067
t_cells <- read.csv('GSE197067_HTSeq_counts.csv', row.names = 1)
rownames(t_cells) <- str_split_fixed(rownames(t_cells), "\\.", 2)[,1]
counts_tpm <- get_tpm(t_cells)
pan_t_metadata <- read.csv('pan_t_metadata.csv', sep = '\t')
pan_t_metadata <- pan_t_metadata[pan_t_metadata$cell_type == 'act',]
counts_tpm <- counts_tpm[pan_t_metadata$cell_type == 'act',]

counts_rictor_tpm <- counts_tpm['ENSG00000164327',]
counts_rbl2_tpm <- counts_tpm['ENSG00000103479',]
counts_usp30_tpm <- counts_tpm['ENSG00000135093',]
counts_hn1_tpm <- counts_tpm['ENSG00000189159',]
counts_osbpl3_tpm <- counts_tpm['ENSG00000070882',]

rictor_expr_d6 <- cbind(pan_t_metadata$sample, pan_t_metadata$cell_type, pan_t_metadata$time, as.vector(counts_rictor_tpm), 'RICTOR')
rbl2_expr_d6 <- cbind(pan_t_metadata$sample, pan_t_metadata$cell_type, pan_t_metadata$time, as.vector(counts_rbl2_tpm), 'RBL2')
usp30_expr_d6 <- cbind(pan_t_metadata$sample, pan_t_metadata$cell_type, pan_t_metadata$time, as.vector(counts_usp30_tpm), 'USP30')
hn1_expr_d6 <- cbind(pan_t_metadata$sample, pan_t_metadata$cell_type, pan_t_metadata$time, as.vector(counts_hn1_tpm), 'HN1')
osbpl3_expr_d6 <- cbind(pan_t_metadata$sample, pan_t_metadata$cell_type, pan_t_metadata$time, as.vector(counts_osbpl3_tpm), 'OSBPL3')

colnames(rictor_expr_d6) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
colnames(rbl2_expr_d6) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
colnames(usp30_expr_d6) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
colnames(hn1_expr_d6) <- c('sample', 'cell_type', 'time', 'expr', 'gene')
colnames(osbpl3_expr_d6) <- c('sample', 'cell_type', 'time', 'expr', 'gene')

rictor_expr_d6 <- as.data.frame(rictor_expr_d6)
rbl2_expr_d6 <- as.data.frame(rbl2_expr_d6)
usp30_expr_d6 <- as.data.frame(usp30_expr_d6)
hn1_expr_d6 <- as.data.frame(hn1_expr_d6)
osbpl3_expr_d6 <- as.data.frame(osbpl3_expr_d6)

rictor_expr_d6[rictor_expr_d6$time == '0h', 'time'] = '0'
rictor_expr_d6[rictor_expr_d6$time == '6h', 'time'] = '6'
rictor_expr_d6[rictor_expr_d6$time == '12h', 'time'] = '12'
rictor_expr_d6[rictor_expr_d6$time == '24h', 'time'] = '24'
rictor_expr_d6[rictor_expr_d6$time == '48h', 'time'] = '48'
rictor_expr_d6[rictor_expr_d6$time == '72h', 'time'] = '72'

rbl2_expr_d6[rbl2_expr_d6$time == '0h', 'time'] = '0'
rbl2_expr_d6[rbl2_expr_d6$time == '6h', 'time'] = '6'
rbl2_expr_d6[rbl2_expr_d6$time == '12h', 'time'] = '12'
rbl2_expr_d6[rbl2_expr_d6$time == '24h', 'time'] = '24'
rbl2_expr_d6[rbl2_expr_d6$time == '48h', 'time'] = '48'
rbl2_expr_d6[rbl2_expr_d6$time == '72h', 'time'] = '72'

usp30_expr_d6[usp30_expr_d6$time == '0h', 'time'] = '0'
usp30_expr_d6[usp30_expr_d6$time == '6h', 'time'] = '6'
usp30_expr_d6[usp30_expr_d6$time == '12h', 'time'] = '12'
usp30_expr_d6[usp30_expr_d6$time == '24h', 'time'] = '24'
usp30_expr_d6[usp30_expr_d6$time == '48h', 'time'] = '48'
usp30_expr_d6[usp30_expr_d6$time == '72h', 'time'] = '72'

hn1_expr_d6[hn1_expr_d6$time == '0h', 'time'] = '0'
hn1_expr_d6[hn1_expr_d6$time == '6h', 'time'] = '6'
hn1_expr_d6[hn1_expr_d6$time == '12h', 'time'] = '12'
hn1_expr_d6[hn1_expr_d6$time == '24h', 'time'] = '24'
hn1_expr_d6[hn1_expr_d6$time == '48h', 'time'] = '48'
hn1_expr_d6[hn1_expr_d6$time == '72h', 'time'] = '72'

osbpl3_expr_d6[osbpl3_expr_d6$time == '0h', 'time'] = '0'
osbpl3_expr_d6[osbpl3_expr_d6$time == '6h', 'time'] = '6'
osbpl3_expr_d6[osbpl3_expr_d6$time == '12h', 'time'] = '12'
osbpl3_expr_d6[osbpl3_expr_d6$time == '24h', 'time'] = '24'
osbpl3_expr_d6[osbpl3_expr_d6$time == '48h', 'time'] = '48'
osbpl3_expr_d6[osbpl3_expr_d6$time == '72h', 'time'] = '72'

rictor_expr_d6$expr <- scale(as.double(rictor_expr_d6$expr))
rbl2_expr_d6$expr <- scale(as.double(rbl2_expr_d6$expr))
usp30_expr_d6$expr <- scale(as.double(usp30_expr_d6$expr))
hn1_expr_d6$expr <- scale(as.double(hn1_expr_d6$expr))
osbpl3_expr_d6$expr <- scale(as.double(osbpl3_expr_d6$expr))

# concat data

rictor_expr <- rbind(rictor_expr_d1, rictor_expr_d2, rictor_expr_d3, rictor_expr_d4, rictor_expr_d5, rictor_expr_d6)#)
rictor_expr$expr <- as.double(rictor_expr$expr)
rictor_expr$time <- factor(as.character(rictor_expr$time), levels = c( "0", "0.5", "1", "2", "4", "6", "8", "12", "24", "48", "72", "96", "120", "144"))

rbl2_expr <- rbind(rbl2_expr_d1, rbl2_expr_d2, rbl2_expr_d3, rbl2_expr_d4, rbl2_expr_d5, rbl2_expr_d6)
rbl2_expr$expr <- as.double(rbl2_expr$expr)
rbl2_expr$time <- factor(as.character(rbl2_expr$time), levels = c( "0", "0.5", "1", "2", "4", "6", "8", "12", "24", "48", "72", "96", "120", "144"))

usp30_expr <- rbind(usp30_expr_d1, usp30_expr_d2, usp30_expr_d3, usp30_expr_d4, usp30_expr_d5, usp30_expr_d6)#)
usp30_expr$expr <- as.double(usp30_expr$expr)
usp30_expr$time <- factor(as.character(usp30_expr$time), levels = c( "0", "0.5", "1", "2", "4", "6", "8", "12", "24", "48", "72", "96", "120", "144"))


hn1_expr <- rbind(hn1_expr_d1, hn1_expr_d2, hn1_expr_d3, hn1_expr_d4, hn1_expr_d5, hn1_expr_d6)
hn1_expr$expr <- as.double(hn1_expr$expr)
hn1_expr$time <- factor(as.character(hn1_expr$time), levels = c( "0", "0.5", "1", "2", "4", "6", "8", "12", "24", "48", "72", "96", "120", "144"))

osbpl3_expr <- rbind(osbpl3_expr_d1, osbpl3_expr_d2, osbpl3_expr_d3, osbpl3_expr_d4, osbpl3_expr_d5, osbpl3_expr_d6)#)
osbpl3_expr$expr <- as.double(osbpl3_expr$expr)
osbpl3_expr$time <- factor(as.character(osbpl3_expr$time), levels = c( "0", "0.5", "1", "2", "4", "6", "8", "12", "24", "48", "72", "96", "120", "144"))

all_expr <- rbind(rictor_expr, rbl2_expr, usp30_expr, hn1_expr, osbpl3_expr)

p <- ggplot(all_expr, aes(x = time, y = expr, fill = gene)) + # , fill = celltype
  facet_wrap(.~gene, ncol = 1) +  # geom_text(label = "test") +
  theme_bw() +
  geom_boxplot() +
  theme(panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  text = element_text(size = 6),
  axis.text.x = element_text(angle = 90, size = 6),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  legend.title = element_blank(), # Title
  plot.title = element_text(size = 6, face = "plain"),
  plot.tag = element_text(size = 8, face = "plain")) +
  geom_smooth(method="loess", aes(group = 1)) +
  labs(x = "Time (hours)" , y = "Expression (TPM)", title = "", tag = "")

ggsave(paste0("timeconc_fantom_intersection.png"), plot = p, units = "mm", width = 170, height = 200)