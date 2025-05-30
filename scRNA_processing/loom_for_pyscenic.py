import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp

#load h5ad
COPD = sc.read("/datf/mazhuo/jupyter_notebook/COPD/processed_data_py/COPD_newdata_0808_celltype.h5ad")
COPD.obs["celltype"] = COPD.obs["celltype"].str.replace(' ', '_')

celltype = ['Mature AT1','AT2','TRB Secretory','PreTB Secretory','Basal','Goblet','Ciliated',
            'Adventitial fibroblast','Alveolar fibroblast','Myofibroblast','Pericyte','ASM','VSM',
            'Mesothelial','Aerocyte','gCap','Arterial','Venous','B','Plasma','CD4+ T','NK','cDC2',
            'Classical monocyte','Non-classical monocyte','Alveolar macrophage',
            'Monocyte-derived macrophage','Interstitial macrophage','Mast','Neutrophil']

#load DEGs
DEG_HCHS = pd.read_csv("/datf/mazhuo/jupyter_notebook/COPD/tables/COPD_DEG_wilcox_HCHS_filter.csv", index_col=0)
DEG_HCHS["celltype"] = DEG_HCHS["celltype"].str.replace(' ', '_')
DEG_HSG12 = pd.read_csv("/datf/mazhuo/jupyter_notebook/COPD/tables/COPD_DEG_wilcox_HSG12_filter.csv", index_col=0)
DEG_HSG12["celltype"] = DEG_HSG12["celltype"].str.replace(' ', '_')
DEG_G1234 = pd.read_csv("/datf/mazhuo/jupyter_notebook/COPD/tables/COPD_DEG_wilcox_G1234_filter.csv", index_col=0)
DEG_G1234["celltype"] = DEG_G1234["celltype"].str.replace(' ', '_')

#filter genes
HCHS = COPD[COPD.obs["Group"].isin(["HC","HS"])].copy()
HSG12 = COPD[COPD.obs["Group"].isin(["HS","Mild COPD"])].copy()
G1234 = COPD[COPD.obs["Group"].isin(["Mild COPD","Severe COPD"])].copy()

#output loom files
for cell in ["Immature_AT1","Mature_AT1","AT2","TRB_Secretory","PreTB_Secretory","Basal","Goblet","Ciliated","PNEC",
                 "Adventitial_fibroblast","Alveolar_fibroblast","Myofibroblast","Fibromyocyte","Pericyte","ASM",
                 "VSM","Mesothelial","Lymphatic","Aerocyte","gCap","Arterial","Venous","B","Plasma","Treg","CD4+_T",
                 "CD8+_T","Proliferating_T","NKT","NK","Basophil","cDC1","cDC2","Migratory_DC","pDC","Classical_monocyte","Non-classical_monocyte",
                 "Alveolar_macrophage","Monocyte-derived_macrophage","Interstitial_macrophage","Mast","Neutrophil"]:
    
    gene = list(set(DEG_HCHS["gene"][DEG_HCHS.celltype == cell].tolist()))
    if len(gene) >= 10:
        cell_ann = HCHS[HCHS.obs["celltype"] == cell].copy()
        cell_ann = cell_ann[:, gene]
        row_attrs = {
            "Gene": np.array(cell_ann.var_names)  
            }
        col_attrs = {
            "CellID": np.array(cell_ann.obs_names)
            }
        if hasattr(cell_ann.X, "toarray"):
            data_matrix = cell_ann.X.toarray()  
        else:
            data_matrix = cell_ann.X
        lp.create("./SCENIC/" + cell + "_DEG_HCHS.loom", data_matrix.transpose(), row_attrs, col_attrs)

for cell in ["Immature_AT1","Mature_AT1","AT2","TRB_Secretory","PreTB_Secretory","Basal","Goblet","Ciliated","PNEC",
                 "Adventitial_fibroblast","Alveolar_fibroblast","Myofibroblast","Fibromyocyte","Pericyte","ASM",
                 "VSM","Mesothelial","Lymphatic","Aerocyte","gCap","Arterial","Venous","B","Plasma","Treg","CD4+_T",
                 "CD8+_T","Proliferating_T","NKT","NK","Basophil","cDC1","cDC2","Migratory_DC","pDC","Classical_monocyte","Non-classical_monocyte",
                 "Alveolar_macrophage","Monocyte-derived_macrophage","Interstitial_macrophage","Mast","Neutrophil"]:
    
    gene = list(set(DEG_HSG12["gene"][DEG_HSG12.celltype == cell].tolist()))
    if len(gene) >= 10:
        cell_ann = HSG12[HSG12.obs["celltype"] == cell].copy()
        cell_ann = cell_ann[:, gene]
        row_attrs = {
            "Gene": np.array(cell_ann.var_names) 
            }
        col_attrs = {
            "CellID": np.array(cell_ann.obs_names)
            }
        if hasattr(cell_ann.X, "toarray"):
            data_matrix = cell_ann.X.toarray()  
        else:
            data_matrix = cell_ann.X
        lp.create("./SCENIC/" + cell + "_DEG_HSG12.loom", data_matrix.transpose(), row_attrs, col_attrs)

for cell in ["Immature_AT1","Mature_AT1","AT2","TRB_Secretory","PreTB_Secretory","Basal","Goblet","Ciliated","PNEC",
                 "Adventitial_fibroblast","Alveolar_fibroblast","Myofibroblast","Fibromyocyte","Pericyte","ASM",
                 "VSM","Mesothelial","Lymphatic","Aerocyte","gCap","Arterial","Venous","B","Plasma","Treg","CD4+_T",
                 "CD8+_T","Proliferating_T","NKT","NK","Basophil","cDC1","cDC2","Migratory_DC","pDC","Classical_monocyte","Non-classical_monocyte",
                 "Alveolar_macrophage","Monocyte-derived_macrophage","Interstitial_macrophage","Mast","Neutrophil"]:
    
    gene = list(set(DEG_G1234["gene"][DEG_G1234.celltype == cell].tolist()))
    if len(gene) >= 10:
        cell_ann = G1234[G1234.obs["celltype"] == cell].copy()
        cell_ann = cell_ann[:, gene]
        row_attrs = {
            "Gene": np.array(cell_ann.var_names) 
            }
        col_attrs = {
            "CellID": np.array(cell_ann.obs_names)
            }
        if hasattr(cell_ann.X, "toarray"):
            data_matrix = cell_ann.X.toarray()  
        else:
            data_matrix = cell_ann.X
        lp.create("./SCENIC/" + cell + "_DEG_G1234.loom", data_matrix.transpose(), row_attrs, col_attrs)
    
