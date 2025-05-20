library(Seurat)
library(SeuratData)
library(ggplot2)
library(RColorBrewer)

data_name <- 'AIDA2'
my_pbmc_s <- readRDS("AIDA_data/ff5a921e-6e6c-49f6-9412-ad9682d23307.rds")
my_pbmc_s <- NormalizeData(my_pbmc_s)
my_pbmc_s <- subset(x = my_pbmc_s, subset = Annotation_Level1 == "T")
DefaultAssay(my_pbmc_s) <- "RNA"
Idents(object = my_pbmc_s) <- my_pbmc_s@meta.data$'Annotation_Level3'
marks <- FindAllMarkers(my_pbmc_s,
  thresh.use = 0,
  min.diff.pct = 0,
  min.pct = 0,
  return.thresh = 1,
  assay = "RNA",
  features = rownames(my_pbmc_s@assays$RNA),
  recorrect_umi = FALSE)

write.table(marks, file = paste0(data_name,"_l3_wilcoxon_marks.csv"), sep="\t", col.names=TRUE, quote = FALSE)
