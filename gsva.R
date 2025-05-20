library(GSVA)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(GO.db)
library(org.Hs.eg.db)

GO <- as.list(GOTERM)

BiocParallel::register(BiocParallel::MulticoreParam(workers = 24, progressbar = TRUE))

# Calculate GSVA scores for each metacell
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
norm_counts <-  as.matrix(obj@assays$RNA@layers$data)
colnames(norm_counts) <- colnames(obj)
rownames(norm_counts) <- rownames(obj)
hs <- org.Hs.eg.db
my.symbols <- rownames(obj)
map_ids <- AnnotationDbi::select(hs, 
       keys = my.symbols,
       columns = c("ENTREZID", "ENSEMBL"),
       keytype = "ENSEMBL")

map_ids <- map_ids[!is.na(map_ids$ENTREZID),]
norm_counts <- norm_counts[map_ids$ENSEMBL, ]
rownames(norm_counts) <- map_ids$ENTREZID
norm_counts <- norm_counts[!duplicated(rownames(norm_counts)),]
goannot <- AnnotationDbi::select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns="GO")
head(goannot)
goannot_bp <- goannot[goannot$ONTOLOGY == 'BP',]
genesbygo <- split(goannot_bp$ENTREZID, goannot_bp$GO)
gsvaPar <- gsvaParam(norm_counts, genesbygo,
                      minSize=5, maxSize=500)
go_res <- gsva(gsvaPar)
saveRDS(go_res, 'gsva_res_all.rds')

# Correlation between CHASERR expression and GSVA scores
expression <- obj@assays$RNA@layers$data
colnames(expression) <- colnames(obj)
rownames(expression) <- rownames(obj)
expression <- t(expression)

x <- as.matrix(expression[,'ENSG00000272888'])
y <- t(as.matrix(go_res))
methods <- c("spearman")
corr_types <- c('rho')
method <- methods[1]
corr_type <- corr_types[1]
goid <- c()
corrs <- c()
pvalues <- c() 
for (i in 1:dim(y)[2])
{
  goid <- c(goid, colnames(y)[i])
  corr <- cor.test(x, y[,i], method = method)
  corrs <- c(corrs, corr$estimate[corr_type]) #
  pvalues <-c(pvalues, corr$p.value)
}

res1 <- data.frame(goid, corrs, pvalues)
res1 <- res1[!is.na(res1$corrs),]
res1$pvalues_adj  <- p.adjust(res1$pvalues, method = 'bonferroni')
res_adj1 <- res1[res1$pvalues_adj < 0.05,]
res_adj1 <- res_adj1[order(abs(res_adj1$corrs), decreasing=TRUE),]
res_adj1$corrs <- round(res_adj1$corrs, digits=4)
anno <- AnnotationDbi::select(GO.db, keys=res_adj1$goid, keytype="GOID", columns=c("TERM","ONTOLOGY") )
rownames(anno) <- anno$GOID
res_adj1$name <- anno[res_adj1$goid, 'TERM']
write.table(res_adj1, file = "gsva_all_corr.csv", sep="\t", col.names=TRUE, quote = FALSE)
