library(ggplot2)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(AnnotationHub)
library(biomaRt)

# A. Polarization data from https://doi.org/10.1038/s41467-020-15543-y
counts <- read.csv('NCOMMS-19-7936188_bulk_RNAseq_raw_counts.txt', sep = '\t')
metadata <- read.csv('NCOMMS-19-7936188_bulk_RNAseq_metadata.txt', sep = '\t')

hub <- AnnotationHub()
query(hub, c("homo sapiens","ensdb"))
ensdb <- hub[["AH116291"]]
require("ensembldb")

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene_lengths =  getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id', 'transcript_length','cds_length'), filters =  'ensembl_gene_id', values = rownames(counts), mart = ensembl, useCache = FALSE)
gene_canonical_transcript =  getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','transcript_is_canonical'), filters = 'ensembl_gene_id', values = rownames(counts), mart = ensembl, useCache = FALSE)
gene_canonical_transcript_subset = gene_canonical_transcript[!is.na(gene_canonical_transcript$transcript_is_canonical),]
gene_lengths = merge(gene_canonical_transcript_subset, gene_lengths, by = c("ensembl_gene_id", "ensembl_transcript_id"))
gene_lengths <- gene_lengths[!is.na(gene_lengths$transcript_length),]
counts <- counts[gene_lengths$ensembl_gene_id,]
x <- counts / gene_lengths$transcript_length
counts_tpm <- t( t(x) * 1e6 / colSums(x) )

counts_chaserr_tpm <- counts_tpm['ENSG00000272888',]
counts_chd2_tpm <- counts_tpm['ENSG00000173575',]

cell <- paste0(metadata$cell_type,   metadata$cytokine_condition,   metadata$stimulation_time)
df_chaserr <- data.frame(log(as.data.frame(counts_chaserr_tpm)), cell, metadata$cell_type,   metadata$cytokine_condition,   metadata$stimulation_time)
colnames(df_chaserr) <- c('value', 'celltype', 'cell_type', 'cytokine_condition', 'stimulation_time')
df_chaserr$gene <- 'CHASERR'
write.table(df_chaserr, file = paste0("chaserr_effectorness.csv"), sep="\t", col.names=TRUE, quote = FALSE)

cell <- paste0(metadata$cell_type,   metadata$cytokine_condition,   metadata$stimulation_time)
df_chd2 <- data.frame(log(as.data.frame(counts_chd2_tpm)), cell, metadata$cell_type,   metadata$cytokine_condition,   metadata$stimulation_time)
colnames(df_chd2) <- c('value', 'celltype', 'cell_type', 'cytokine_condition', 'stimulation_time')
df_chd2$gene <- 'CHD2'
write.table(df_chd2, file = paste0("chd2_effectorness.csv"), sep="\t", col.names=TRUE, quote = FALSE)

# plot
effectorness_chaserr = read.csv('chaserr_effectorness.csv', sep = '\t')
effectorness_chaserr$stimulation_time <- factor(effectorness_chaserr$stimulation_time, levels = c('16h', '5d'))
effectorness_chaserr[effectorness_chaserr$cell_type == 'CD4_Naive', 'cell_type'] = 'CD4+ naive'
effectorness_chaserr[effectorness_chaserr$cell_type == 'CD4_Memory', 'cell_type'] = 'CD4+ activated'
effectorness_chaserr$cell_type <- factor(effectorness_chaserr$cell_type, levels = c('CD4+ naive', 'CD4+ activated'))
effectorness_chaserr$cytokine_condition <- factor(effectorness_chaserr$cytokine_condition, levels = c("Resting", "iTreg", "Th17", "IFNB", "Th0", "Th1", "Th2"))

p1 <- ggplot(effectorness_chaserr, aes(x = cytokine_condition, y = value, fill=cytokine_condition))+ 
  facet_grid(~ cell_type + stimulation_time) +	
  theme_bw() + 
  geom_boxplot(fatten = 0.8, outlier.shape = NA, alpha = 0.4) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 6, face = "plain"),
        plot.tag = element_text(size = 8, face = "plain")) +
    labs(y = "Expression (TPM)", title = "CHASERR", tag = "A")

effectorness_chd2 = read.csv('chd2_effectorness.csv', sep = '\t')
effectorness_chd2$stimulation_time <- factor(effectorness_chd2$stimulation_time, levels = c('16h', '5d'))

effectorness_chd2[effectorness_chd2$cell_type == 'CD4_Naive', 'cell_type'] = 'CD4+ naive'
effectorness_chd2[effectorness_chd2$cell_type == 'CD4_Memory', 'cell_type'] = 'CD4+ activated'
effectorness_chd2$cell_type <- factor(effectorness_chd2$cell_type, levels = c('CD4+ naive', 'CD4+ activated'))
effectorness_chd2$cytokine_condition <- factor(effectorness_chd2$cytokine_condition, levels = c("Resting", "iTreg", "Th17", "IFNB", "Th0", "Th1", "Th2"))

p2 <- ggplot(effectorness_chd2, aes(x = cytokine_condition, y = value, fill=cytokine_condition))+ 
  facet_grid(~ cell_type + stimulation_time) +	
  theme_bw() + 
  geom_boxplot(fatten = 0.8, outlier.shape = NA, alpha = 0.4) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 6, face = "plain"),
        plot.tag = element_text(size = 8, face = "plain")) +
    labs(y = "Expression (TPM)", title = "CHD2", tag = "")

p <-ggarrange(p1, p2, ncol=1, nrow=2, common.legend = TRUE, legend="right")


p3 <- readRDS('time_activation_main.rds')
g <- arrangeGrob(p, p3, ncol = 1, nrow = 2, layout_matrix= rbind(c(1), c(2)))
ggsave(paste0("cd4_polarization.png"), plot = g, units = "mm", width = 170, height = 200)
