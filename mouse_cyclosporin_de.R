library(Seurat)
library(RCAv2)
library(RColorBrewer)
library(ggplot2)
library(stringr)

counts_t_s <- readRDS("mouse2_T.rds")
counts_t_s_ctl <- subset(counts_t_s, sample %in% c("CTL1", "CTL2", "CTL3", "CTL4"))
counts_t_s_eau <- subset(counts_t_s, sample %in% c("EAU1", "EAU2", "EAU3", "EAU4"))
counts_t_s_csa <- subset(counts_t_s, sample %in% c("CSA1", "CSA2", "CSA3", "CSA4"))

DefaultAssay(counts_t_s_ctl) <- "RNA"
Idents(object = counts_t_s_ctl) <- counts_t_s_ctl@meta.data$celltype_l2
marks_ctl <- FindAllMarkers(counts_t_s_ctl, return.thresh = 1, assay = "RNA", features = rownames(counts_t_s_ctl@assays$RNA), recorrect_umi = FALSE)
marks_ctl$state <- 'CTL'

DefaultAssay(counts_t_s_eau) <- "RNA"
Idents(object = counts_t_s_eau) <- counts_t_s_eau@meta.data$celltype_l2
marks_eau <- FindAllMarkers(counts_t_s_eau, return.thresh = 1, assay = "RNA", features = rownames(counts_t_s_eau@assays$RNA), recorrect_umi = FALSE)
marks_eau$state <- 'EAU'

DefaultAssay(counts_t_s_csa) <- "RNA"
Idents(object = counts_t_s_csa) <- counts_t_s_csa@meta.data$celltype_l2
marks_csa <- FindAllMarkers(counts_t_s_csa, return.thresh = 1, assay = "RNA", features = rownames(counts_t_s_csa@assays$RNA), recorrect_umi = FALSE)
marks_csa$state <- 'CSA'

marks <- rbind(marks_ctl, marks_eau, marks_csa)
write.table(marks, file = paste0("cycl_wilcoxon_marks.csv"), sep="\t", col.names=TRUE, quote = FALSE)

group_cell <- counts_t_s@meta.data$sample
group_cell[group_cell %in% c('EAU1','EAU2','EAU3','EAU4')] <- 'EAU'
group_cell[group_cell %in% c('CTL1','CTL2','CTL3','CTL4')] <- 'CTL'
group_cell[group_cell %in% c('CSA1','CSA2','CSA3','CSA4')] <- 'CSA'

counts_t_s@meta.data$de_group <- group_cell
Idents(object = counts_t_s) <- counts_t_s@meta.data$de_group
celltypes <- c("CD4 naive", "Treg", "CD8 naive", "CD4 Th17", "CD4 Tfh", "CD4 Th1", "Pro T", "CD8 CTL")
init <- FALSE
for (ct in celltypes)
{
  counts_s_cluster <- subset(x = counts_t_s, subset = celltype_l2 == ct)
  marks1 <- FindMarkers(counts_s_cluster, return.thresh = 1, assay = "RNA", ident.1 = "EAU",  ident.2 = "CTL" , features = rownames(counts_s_cluster@assays$RNA), recorrect_umi = FALSE) #, test.use = 'MAST'
  marks1$celltype <- ct
  marks1$comp <- 'CTL_EAU'
  marks1$gene <- rownames(marks1)
  if(!init)
  {
    init <- TRUE
    marks <- marks1
  }
  else
  {
    marks <- rbind(marks, marks1)
  }
  
  marks1 <- FindMarkers(counts_s_cluster, return.thresh = 1, assay = "RNA", ident.1 = "CSA" ,  ident.2 = "EAU", features = rownames(counts_s_cluster@assays$RNA), recorrect_umi = FALSE) #, test.use = 'MAST'
  marks1$celltype <- ct
  marks1$comp <- 'EAU_CSA'
  marks1$gene <- rownames(marks1)
  marks <- rbind(marks, marks1)
}

write.table(marks, paste0("marks_CTL_EAU_CSA.csv"), row.names = F, quote = F, sep="\t")

