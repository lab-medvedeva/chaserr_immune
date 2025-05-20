library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library("org.Hs.eg.db")
library(BiocParallel)
library(clusterProfiler)
library(AnnotationHub)
library(biomaRt)
require("ensembldb")

enableWGCNAThreads(nThreads = 16)

hub <- AnnotationHub()
query(hub, c("homo sapiens","ensdb"))
ensdb <- hub[["AH116291"]]
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")


hs <- org.Hs.eg.db

my_pbmc_s <- readRDS('aida_metacells.rds')
obj <- my_pbmc_s@misc$tutorial$wgcna_metacell_obj
obj <- subset(x = obj, subset = cell_type !=  'CD4-positive, alpha-beta T cell')
obj <- subset(x = obj, subset = cell_type !=  'CD8-positive, alpha-beta T cell')
obj[["RNA"]]$data <- obj[["RNA"]]$counts
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:30)

universe <- rownames(obj)

for(ct in c('naive thymus-derived CD4-positive, alpha-beta T cell', 
            'naive thymus-derived CD8-positive, alpha-beta T cell',
            'central memory CD4-positive, alpha-beta T cell',
            'effector memory CD4-positive, alpha-beta T cell',
            'CD4-positive, alpha-beta cytotoxic T cell',
            'CD8-positive, alpha-beta cytotoxic T cell',
            "CD8-positive, alpha-beta memory T cell",
            'regulatory T cell'))
{
  obj_type <- subset(x = obj, subset = cell_type == ct)
  expression <- obj_type@assays$RNA@layers$data
  colnames(expression) <- colnames(obj_type)
  rownames(expression) <- rownames(obj_type)
  expression <- t(expression)

  # Spearman correlation
  x <- as.matrix(expression[,'ENSG00000272888'])
  y <- as.matrix(expression)
  methods <- c("spearman")
  corr_types <- c('rho')
  method <- methods[1]
  corr_type <- corr_types[1]
  genes <- c()
  corrs <- c()
  pvalues <- c() 
  for (i in 1:dim(y)[2])
  {
    genes <- c(genes, colnames(expression)[i])
    corr <- cor.test(x, y[,i], method = method)
    corrs <- c(corrs, corr$estimate[corr_type]) #
    pvalues <- c(pvalues, corr$p.value)
  }

  res1 <- data.frame(genes, corrs, pvalues)
  res1 <- res1[!is.na(res1$corrs),]
  res1$pvalues_adj  <- p.adjust(res1$pvalues, method = 'bonferroni')
  res_adj1 <- res1[res1$pvalues_adj < 0.05,]
  res_adj1 <- res_adj1[order(abs(res_adj1$corrs), decreasing=TRUE),]
  res_adj_cd4 <- res_adj1[abs(res_adj1$corrs) > 0.5,]

  # # Convert to gene symbols
  my.symbols <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), 
                            values= res_adj_cd4$genes, 
                            filters = "ensembl_gene_id", 
                            mart = mart, 
                            uniqueRows = TRUE)
  my.symbols <- my.symbols[!duplicated(my.symbols$hgnc_symbol), ]
  my.symbols <- my.symbols[!is.na(my.symbols$hgnc_symbol), ]
  my.symbols <- my.symbols[my.symbols$hgnc_symbol!="", ]

  rownames(my.symbols) <- my.symbols$ensembl_gene_id
  res_adj_cd4$Symbol <- my.symbols[res_adj_cd4$gene, 'hgnc_symbol']
  write.table(res_adj_cd4, file = paste0(ct, "_corr.csv"), sep="\t", col.names=TRUE, row.names=FALSE, quote = FALSE)
  
  # GO:BP and GSEA
  res_adj_cd4 <- read.csv(file = paste0(ct, "_corr.csv"), sep="\t")
  res_adj_cd4_pos <- dplyr::filter(res_adj_cd4, res_adj_cd4$corrs > 0 )
  res_adj_cd4_neg <- dplyr::filter(res_adj_cd4, res_adj_cd4$corrs < 0 )
  ego_BP_pos <- enrichGO(gene      = res_adj_cd4_pos$genes,
                     OrgDb         = org.Hs.eg.db,
                     universe      = universe,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05)

  write.table(as.data.frame(ego_BP_pos), file = paste0(ct, "_go_bp_pos.csv"), sep="\t", col.names=TRUE, quote = FALSE)
  saveRDS(ego_BP_pos,  paste0(ct, "_go_bp_pos.rds"))
  ego_BP_neg <- enrichGO(gene      = res_adj_cd4_neg$genes,
                     OrgDb         = org.Hs.eg.db,
                     universe      = universe,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05)

  write.table(as.data.frame(ego_BP_neg), file = paste0(ct, "_go_bp_neg.csv"), sep="\t", col.names=TRUE, quote = FALSE)
  saveRDS(ego_BP_neg,  paste0(ct, "_go_bp_neg.rds"))
  df <- res_adj_cd4[order(res_adj_cd4$corrs, decreasing=TRUE),]
  geneList <- df$corrs
  names(geneList) <- df$genes
  gsea <- gseGO(gene = geneList,
             OrgDb         = org.Hs.eg.db, 
             keyType       = 'ENSEMBL',
             ont           = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff   = 0.1)
  
  write.table(as.data.frame(gsea), file = paste0(ct, "_gsea.csv"), sep="\t", col.names=TRUE, quote = FALSE)
  saveRDS(gsea,  paste0(ct, "_gsea.rds"))
}

