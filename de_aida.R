library(Seurat)
library(SeuratData)
library(ggplot2)
library(RColorBrewer)
library(splitstackshape)
library(stringr)
library(DESeq2)

data_name <- 'AIDA2'
my_pbmc_s <- readRDS("AIDA_data/ff5a921e-6e6c-49f6-9412-ad9682d23307.rds")
my_pbmc_s <- NormalizeData(my_pbmc_s)

# AIDA, compare T with other cell types from L1 annotation
pseudo_my_pbmc_s <- AggregateExpression(my_pbmc_s, assays = "RNA", return.seurat = T, group.by = c("donor_id", "Annotation_Level1"))
pseudo_my_pbmc_s@meta.data$cell_type <-   str_split_fixed(rownames(pseudo_my_pbmc_s@meta.data), "_", 2)[,2]
Idents(pseudo_my_pbmc_s) <- "cell_type"
res = list()
for (cell_type in unique(pseudo_my_pbmc_s@meta.data$'cell_type'))
{
	if(cell_type != 'T')
	{
		marks <- FindMarkers(pseudo_my_pbmc_s,
		  ident.1 = cell_type,
		  ident.2 = 'T',
		  logfc.threshold = 0,
		  thresh.use = 0,
		  min.pct = 0.01,
		  return.thresh = 1,
		  assay = "RNA",
		  recorrect_umi = FALSE)
		marks$gene <- rownames(marks)
		marks$cell_type <- cell_type
		res[[cell_type]] <- marks[c('ENSG00000272888', 'ENSG00000173575'),]
	}
}

df <- do.call("rbind", res)
write.table(df, file = paste0(data_name, "_l1_wilcoxon_T.csv"), sep="\t", col.names=TRUE, row.names=TRUE, quote = FALSE)

# AIDA, Differential expression analysis for T cell subtypes, pseudobulk data L3
my_pbmc_s <- subset(x = my_pbmc_s, subset = Annotation_Level1 == "T")
DefaultAssay(my_pbmc_s) <- "RNA"
Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$'Annotation_Level3'
pseudo_my_pbmc_s <- AggregateExpression(my_pbmc_s, assays = "RNA", return.seurat = T, group.by = c("donor_id", "Annotation_Level3"))
pseudo_my_pbmc_s@meta.data$cell_type <-   str_split_fixed(rownames(pseudo_my_pbmc_s@meta.data), "_", 2)[,2]
Idents(pseudo_my_pbmc_s) <- "cell_type"
marks <- FindAllMarkers(pseudo_my_pbmc_s,
  thresh.use = 0,
  logfc.threshold = 0,
  min.pct = 0.01,
  return.thresh = 1,
  assay = "RNA",
  recorrect_umi = FALSE)

write.table(marks, file = paste0(data_name,"_l3_wilcoxon_marks_pseudobulk.csv"), sep="\t", col.names=TRUE, quote = FALSE)

# AIDA, Differential expression analysis for T cell subtypes, single-cell data L3
marks <- FindAllMarkers(my_pbmc_s,
  thresh.use = 0,
  logfc.threshold = 0,
  min.pct = 0.01,
  return.thresh = 1,
  assay = "RNA",
  recorrect_umi = FALSE)

write.table(marks, file = paste0(data_name,"_l3_wilcoxon_marks.csv"), sep="\t", col.names=TRUE, quote = FALSE)

