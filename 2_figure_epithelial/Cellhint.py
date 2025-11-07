import scanpy as sc
import anndata
import cellhint 
import seaborn as sns
import pandas as pd
import numpy as np

all = sc.read("./processed_data_py/COPD_all_cellhint.h5ad")
epi = all[all.obs["cellclass"] == "Epithelial"].copy()

sc.pp.highly_variable_genes(epi, n_top_genes=5000, batch_key = 'dataset', flavor = "seurat")
epi = epi[:, epi.var.highly_variable] 
sc.pp.scale(epi)
sc.pp.regress_out(epi, ['dataset'])
sc.tl.pca(epi)
alignment = cellhint.harmonize(epi, 'dataset', 'celltype')
alignment.best_align(dataset_order = ["Our", "Murthy et al., 2022",
                                      "Guo et al., 2022","Adams et al., 2020",
                                      "Travaglini et al. 2020"])
cellhint.treeplot(alignment, figsize=(13,8), save = "./figures/Cellhint_tree_epi.pdf")
