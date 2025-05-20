library(Seurat)
library(RCAv2)
library(RColorBrewer)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(harmony)

# Merge and intergrate samples
samples <- c("EAU1", "EAU2", "EAU3", "EAU4", "CSA1", "CSA2", "CSA3", "CSA4", "CTL1", "CTL2", "CTL3", "CTL4")
sample_rds <- list()
for(sample in samples)
{
  rca_object <- createRCAObjectFrom10X(paste0(sample, "/filtered_feature_bc_matrix"))

  rca_object <- dataFilter(rca_object,
                       nGene.thresholds = c(200, 2500),
                       nUMI.thresholds = c(500, 10000),
                       percent.mito.thresholds = c(0, 0.15), 
                       min.cell.exp = 0.01*ncol(rca_object$raw.data), 
                       plot = FALSE)

  gene.row <- rownames(rca_object$raw.data)
  gene.row.good <- unique(gene.row)
  raw.data.old <- rca_object$raw.data
  rca_object$raw.data <- rca_object$raw.data[gene.row.good,]
  seurat_object <- CreateSeuratObject(counts = rca_object$raw.data, project = "csa", min.cells = 0, min.features = 0)
  colnames(seurat_object) <- paste0(sample, "_", colnames(seurat_object))
  sample_rds[[sample]] <- seurat_object
}

merged_combined <- merge(sample_rds[["EAU1"]], sample_rds[["EAU2"]], merge.data = TRUE)
for (i in 3:length(samples))
{
  merged_combined <- merge(merged_combined, sample_rds[[samples[i]]], merge.data = TRUE)
}

counts_s <- merged_combined
counts_s <- NormalizeData(counts_s)
counts_s <- FindVariableFeatures(counts_s)
counts_s <- ScaleData(counts_s)
counts_s <- RunPCA(counts_s)
counts_s <- FindNeighbors(counts_s, dims = 1:30, reduction = "pca")
counts_s <- FindClusters(counts_s)
counts_s <- RunUMAP(counts_s, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
counts_s <- IntegrateLayers(object = counts_s, method = HarmonyIntegration, orig.reduction = "pca",  new.reduction = 'harmony', verbose = FALSE)
counts_s <- FindNeighbors(counts_s, reduction = "harmony", dims = 1:30)
counts_s <- FindClusters(counts_s, resolution = 0.5, cluster.name = "harmony_clusters")
counts_s <- RunUMAP(counts_s, dims = 1:30, reduction = "harmony", reduction.name = "harmony_umap")
counts_s <- JoinLayers(counts_s)

p <- DimPlot(counts_s, reduction = "harmony_umap", group.by = "harmony_clusters", label = TRUE, label.size = 6, repel = TRUE, raster=FALSE) + theme(legend.position = "bottom", legend.text=element_text(size=12))
ggsave(paste0("mouse_harmony.png"), plot = p, width = 15, height = 15)

counts_s@meta.data$sample <- str_split_fixed(colnames(counts_s), "_", 2)[,1]

# Annotate T cells
DotPlot(
  counts_s,
  features = c('Cd3e', 'Cd3d', 'Cd3g'),
  group.by = "harmony_clusters")

DotPlot(
  counts_s,
  features = c('Cd79a', 'Cd79b', 'Ms4a1'),
  group.by = "harmony_clusters")

DotPlot(
  counts_s,
  features = c('Hbb-a1', 'Hbb-bs'),
  group.by = "harmony_clusters")

counts_s@meta.data$celltype <- NA
counts_s@meta.data$celltype[counts_s@meta.data$harmony_clusters == 1] = "T"
counts_s@meta.data$celltype[counts_s@meta.data$harmony_clusters == 2] = "T"
counts_s@meta.data$celltype[counts_s@meta.data$harmony_clusters == 5] = "T"
counts_s@meta.data$celltype[counts_s@meta.data$harmony_clusters == 8] = "T"
counts_s@meta.data$celltype[counts_s@meta.data$harmony_clusters == 9] = "T"
counts_s@meta.data$celltype[counts_s@meta.data$harmony_clusters == 13] = "T"

p <- DimPlot(counts_s, reduction = "harmony_umap", group.by = "celltype", label = TRUE, label.size = 6, repel = TRUE, raster=FALSE) + theme(legend.position = "bottom", legend.text=element_text(size=12))
ggsave(paste0("mouse2_harmony.png"), plot = p, width = 15, height = 15)
saveRDS(counts_s, "mouse2.rds")

# counts_s <- readRDS("mouse2.rds")

# Recluster T cells
counts_t_s <- subset(counts_s, celltype %in% c("T"))
counts_t_s <- FindClusters(counts_t_s, resolution = 1, cluster.name = "harmony_clusters2")
counts_t_s <- RunUMAP(counts_t_s, dims = 1:30, reduction = "harmony", reduction.name = "harmony_umap2")

p1 <- DimPlot(counts_t_s, reduction = "harmony_umap2", group.by = "harmony_clusters", label = TRUE, label.size = 6, repel = TRUE, raster=FALSE) + theme(legend.position = "bottom", legend.text=element_text(size=12))
p2 <- DimPlot(counts_t_s, reduction = "harmony_umap2", group.by = "harmony_clusters2", label = TRUE, label.size = 6, repel = TRUE, raster=FALSE) + theme(legend.position = "bottom", legend.text=element_text(size=12))
p1|p2

group_cell <- counts_t_s@meta.data$sample
group_cell[group_cell %in% c('EAU1','EAU2','EAU3','EAU4')] <- 'EAU'
group_cell[group_cell %in% c('CTL1','CTL2','CTL3','CTL4')] <- 'CTL'
group_cell[group_cell %in% c('CSA1','CSA2','CSA3','CSA4')] <- 'CSA'

counts_t_s@meta.data$de_group <- group_cell

p1 <- DimPlot(counts_t_s, reduction = "harmony_umap2", group.by = "harmony_clusters2", label = TRUE, label.size = 6, repel = TRUE, raster=FALSE) + theme(legend.position = "bottom", legend.text=element_text(size=12))
p2 <- DimPlot(counts_t_s, reduction = "harmony_umap2", group.by = "harmony_clusters2", label = TRUE, label.size = 6, repel = TRUE, raster=FALSE) + theme(legend.position = "bottom", legend.text=element_text(size=12))
p1|p2

# Annotate T cells
DotPlot(
  counts_t_s,
  features = c('Foxp3','Gzmk','Cd8a','Igfbp4','Gzmb','Ctla2a','Mki67','Stmn1','Sostdc1','Cxcr5','Cd40lg','Cd4','Ccr7','Cxcr6', 'Ifngr1','Il17a', 'Ccr2','Rorc'),
  group.by = "harmony_clusters2")

counts_t_s@meta.data$celltype_l2 <- NA
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 7] = "Treg"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 12] = "Treg"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 1] = "CD8 naive"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 4] = "CD8 naive"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 2] = "CD8 naive"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 5] = "CD8 CTL"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 15] = "CD8 CTL"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 0] = "CD4 naive"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 3] = "CD4 naive"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 6] = "CD4 naive"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 16] = "CD4 naive"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 9] = "CD4 Tfh"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 14] = "CD4 Th17"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 8] = "CD4 Th1"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 11] = "Pro T"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 18] = "CD8 CTL"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 10] = "Other T"
counts_t_s@meta.data$celltype_l2[counts_t_s@meta.data$harmony_clusters2 == 13] = "Other T"

p <- DimPlot(counts_t_s, reduction = "harmony_umap2", group.by = "celltype_l2", label = TRUE, label.size = 6, repel = TRUE, raster=FALSE) + theme(legend.position = "bottom", legend.text=element_text(size=12))
p

saveRDS(counts_t_s, "mouse2_T.rds")
