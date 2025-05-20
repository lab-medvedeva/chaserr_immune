library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library("org.Hs.eg.db")

enableWGCNAThreads(nThreads = 16)

my_pbmc <- readRDS("AIDA_data/ff5a921e-6e6c-49f6-9412-ad9682d23307.rds")
my_pbmc <- NormalizeData(my_pbmc)
my_pbmc_s <- subset(x = my_pbmc, subset = cell_type %in% c('central memory CD4-positive, alpha-beta T cell',
                                'naive thymus-derived CD4-positive, alpha-beta T cell',
                                'effector memory CD4-positive, alpha-beta T cell',
                                'CD4-positive, alpha-beta T cell',
                                'CD4-positive, alpha-beta cytotoxic T cell',
                                'regulatory T cell',
                                'CD8-positive, alpha-beta memory T cell',
                                'naive thymus-derived CD8-positive, alpha-beta T cell',
                                'CD8-positive, alpha-beta cytotoxic T cell',
                                'CD8-positive, alpha-beta T cell'))

tmp_object <- subset(my_pbmc_s, features =  c('ENSG00000272888', 'ENSG00000173575'))
tmp_object[["RNA"]]$data <- as(object = tmp_object[["RNA"]]$data, Class = "dgCMatrix")
Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$library_uuid
m <- AverageExpression(
  tmp_object,
  group.by = "library_uuid",
  features = c('ENSG00000272888', 'ENSG00000173575'),
  assays = 'RNA',
  layer = "data")

df <- as.data.frame(m$RNA)
libraries <- colnames(m$RNA)
chaserr_mean <- as.double(df['ENSG00000272888',])
chd2_mean <- as.double(df['ENSG00000272888',])
five_quantile_chaserr <- quantile(chaserr_mean , probs = 0.1, na.rm = FALSE)[[1]]
five_quantile_chd2 <- quantile(chd2_mean, probs = 0.1, na.rm = FALSE)[[1]]
libraries_filtered <- libraries[(chaserr_mean > five_quantile_chaserr) & (chd2_mean > five_quantile_chd2)]
my_pbmc_s <- subset(x = my_pbmc_s, subset = library_uuid %in% libraries_filtered)
Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$cell_type
my_pbmc_s[["RNA"]] <- split(my_pbmc_s[["RNA"]], f = my_pbmc_s$library_uuid)
my_pbmc_s <- NormalizeData(my_pbmc_s)
my_pbmc_s <- FindVariableFeatures(my_pbmc_s)
my_pbmc_s <- ScaleData(my_pbmc_s)
my_pbmc_s <- RunPCA(my_pbmc_s)
my_pbmc_s <- IntegrateLayers(object = my_pbmc_s, method = HarmonyIntegration, orig.reduction = "pca",  new.reduction = 'harmony', verbose = FALSE)
my_pbmc_s <- FindNeighbors(my_pbmc_s, reduction = "harmony", dims = 1:30)
my_pbmc_s <- FindClusters(my_pbmc_s, resolution = 0.5, cluster.name = "harmony_clusters")
my_pbmc_s <- RunUMAP(my_pbmc_s, dims = 1:30, reduction = "harmony", reduction.name = "harmony_umap")
my_pbmc_s <- JoinLayers(my_pbmc_s)
saveRDS(my_pbmc_s , paste0('AIDA_reintegrated.rds'))

my_pbmc_s <- SetupForWGCNA(
  my_pbmc_s,
  features = rownames(my_pbmc_s),
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "tutorial"
)

my_pbmc_s <- MetacellsByGroups(
  seurat_obj = my_pbmc_s,
  reduction = "harmony",
  slot = "data",
  layer = "data",
  mode = "average",
  group.by = c('cell_type'),
  k = 75,
  max_shared = 10,
  ident.group = 'cell_type'
)
saveRDS(my_pbmc_s, paste0('aida_metacells.rds'))
