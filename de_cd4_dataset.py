import numpy as np
import scanpy as sc
import pandas as pd
import scipy

# CD4 T data https://doi.org/10.1016/j.xgen.2023.100473
data_name = 'Validation'
adata = sc.read_h5ad("res.cxg.h5ad")
sc.tl.rank_genes_groups(adata, groupby = 'cluster_L1', key='wilcoxon')
marker_df = sc.get.rank_genes_groups_df(adata, group = None)
marker_df.to_csv(f"{data_name}_l1_wilcoxon_marks.csv", sep="\t")

sc.tl.rank_genes_groups(adata, groupby = 'cluster_L2', key='wilcoxon')
marker_df = sc.get.rank_genes_groups_df(adata, group = None)
marker_df.to_csv(f"{data_name}_l2_wilcoxon_marks.csv", sep="\t")
