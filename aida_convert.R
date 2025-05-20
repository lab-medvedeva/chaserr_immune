library(Seurat)
library(BPCells)

data <- open_matrix_anndata_hdf5("ff5a921e-6e6c-49f6-9412-ad9682d23307.h5ad")

 write_matrix_dir(
   mat = data,
   dir = "ff5a921e-6e6c-49f6-9412-ad9682d23307"
 )

mat <- open_matrix_dir(dir = "ff5a921e-6e6c-49f6-9412-ad9682d23307")
metadata <- read.csv("ff5a921e-6e6c-49f6-9412-ad9682d23307.metadata", sep = '\t', row.names = 1) 
merged.object <- CreateSeuratObject(counts = mat, meta.data = metadata)
UMAP_coordinates_mat <- read.csv("ff5a921e-6e6c-49f6-9412-ad9682d23307.umap", sep = '\t', header = F)#, row.names = 1) 
UMAP_coordinates_mat <- as.matrix(UMAP_coordinates_mat)
rownames(UMAP_coordinates_mat) <- rownames(metadata)
merged.object[['UMAP']] <- CreateDimReducObject(embeddings = UMAP_coordinates_mat, key = "UMAP_", global = T, assay = "RNA")
saveRDS(merged.object, "AIDA_data/ff5a921e-6e6c-49f6-9412-ad9682d23307.rds")