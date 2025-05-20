library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(stringr)
library(clusterProfiler)
library("org.Mm.eg.db") 

counts_t_s <- readRDS("mouse2_T.rds")

counts_t_s@meta.data$cells <- rownames(counts_t_s@meta.data)
counts_t_s <- subset(x = counts_t_s, cells %in% counts_t_s@meta.data[!is.na(counts_t_s@meta.data$celltype_l2),'cells'])
subset_cells <- WhichCells(counts_t_s, expression = Chaserr > 0, accept.low = 0.1, use.raw=T)
counts_t_s_filtered <- subset(x = counts_t_s, cells %in% subset_cells)
norm_data <- counts_t_s_filtered@assays$RNA@layers$data
rownames(norm_data) <- rownames(counts_t_s_filtered@assays$RNA) 
norm_data <- norm_data[!startsWith(rownames(norm_data), "Rps"),]
norm_data <- norm_data[!startsWith(rownames(norm_data), "Rpl"),]
norm_data <- norm_data[!startsWith(rownames(norm_data), "mt-"),]
x <- log10(10^as.matrix(as.data.frame(norm_data['Chaserr',])) - 1 + 1.01)
y <- log10(10^t(as.matrix(norm_data)) - 1 + 1.01)

methods <- c("spearman")
corr_types <- c('rho')
method <- methods[1]
corr_type <- corr_types[1]
genes <- c()
corrs <- c()
pvalues <- c() 
for (i in 1:dim(y)[2])
{
genes <- c(genes, rownames(norm_data)[i])
corr <- cor.test(x, y[,i], method = method)
corrs <- c(corrs, corr$estimate[corr_type]) #
pvalues <-c(pvalues, corr$p.value)
}

res1 <- data.frame(genes, corrs, pvalues)
res1 <- res1[!is.na(res1$corrs),]
res1$pvalues_adj  <- p.adjust(res1$pvalues, method = 'bonferroni')
res_adj1 <- res1[res1$pvalues_adj < 0.05,]
res_adj1 <- res_adj1[order(res_adj1$corrs, decreasing=TRUE),]
write.table(res_adj1, file = paste0("mouse_cyclosporin_corr.csv"), sep="\t", row.names = FALSE, col.names=TRUE, quote = FALSE)

res_adj1_pos <- res_adj1[res_adj1$corrs > 0,]
ego_BP_pos <- enrichGO(gene      = res_adj1_pos$genes,
                   OrgDb         = org.Mm.eg.db,
                   universe      = rownames(counts_t_s_filtered),
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05)

write.table(as.data.frame(ego_BP_pos), file = paste0("mouse_cyclosporin_go_bp_pos.csv"), sep="\t", row.names = FALSE, col.names=TRUE, quote = FALSE)

res_adj1_neg <- res_adj1[res_adj1$corrs < 0,]
ego_BP_neg <- enrichGO(gene      = res_adj1_neg$genes,
                   OrgDb         = org.Mm.eg.db,
                   universe      = rownames(counts_t_s_filtered),
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05)

write.table(as.data.frame(ego_BP_neg), file = paste0("mouse_cyclosporin_go_bp_neg.csv"), sep="\t", row.names = FALSE, col.names=TRUE, quote = FALSE)


