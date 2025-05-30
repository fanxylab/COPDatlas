library(tidyverse)
library(sscVis)
library(data.table)
library(Startrac)
library(ComplexHeatmap)
library(colorRamp2)

#load cell metadata
meta = read.csv("/datf/mazhuo/jupyter_notebook/COPD/tables/COPD_meta_0808.csv")

#Calculate Ro/e value by each cell class
meta$celltype = as.character(meta$celltype)
meta$sample = as.factor(meta$sample)
epi_meta = meta[which(meta$celltype %in% c("Immature AT1","Mature AT1","AT2","TRB Secretory","PreTB Secretory","Basal","Goblet","Ciliated","PNEC")),]
epi_meta$celltype = factor(epi_meta$celltype, levels = c("Immature AT1","Mature AT1","AT2","TRB Secretory","PreTB Secretory","Basal","Goblet","Ciliated","PNEC"))
startrac.dist.epi = calTissueDist(epi_meta,byPatient = F,
                              colname.cluster = "celltype",
                              colname.patient = "sample",
                              colname.tissue = "new_group")

mes_meta = meta[which(meta$celltype %in% c("Adventitial fibroblast","Alveolar fibroblast","Myofibroblast","Fibromyocyte","Pericyte","ASM","VSM")),]
mes_meta$celltype = factor(mes_meta$celltype, levels = c("Adventitial fibroblast","Alveolar fibroblast","Myofibroblast","Fibromyocyte","Pericyte","ASM","VSM"))
startrac.dist.mes = calTissueDist(mes_meta,byPatient = F,
                              colname.cluster = "celltype",
                              colname.patient = "sample",
                              colname.tissue = "new_group")

endo_meta = meta[which(meta$cellclass == "Endothelial"),]
endo_meta$celltype = factor(endo_meta$celltype, levels = c("Lymphatic","Aerocyte","gCap","Arterial","Venous"))
startrac.dist.endo = calTissueDist(endo_meta,byPatient = F,
                              colname.cluster = "celltype",
                              colname.patient = "sample",
                              colname.tissue = "new_group")

immune_meta = meta[which(meta$cellclass == "Immune"),]
immune_meta$celltype = factor(immune_meta$celltype, levels = c("B","Plasma","Treg","CD4+ T","CD8+ T","Proliferating T","NKT","NK",
                                                               "cDC1","cDC2","Migratory DC","pDC","Classical monocyte",
                                                               "Non-classical monocyte","Alveolar macrophage","Monocyte-derived macrophage",
                                                               "Interstitial macrophage","Mast","Neutrophil","Basophil"))
startrac.dist.immune = calTissueDist(immune_meta,byPatient = F,
                              colname.cluster = "celltype",
                              colname.patient = "sample",
                              colname.tissue = "new_group")

#Combine Ro/e values from each cell class
startrac.dist = rbind(startrac.dist.epi, startrac.dist.mes, startrac.dist.endo, startrac.dist.immune)
startrac.dist = t(round(startrac.dist, 2))

split = data.frame(celltype = colnames(startrac.dist), class = c(rep("Epithelial", 9),
                                                                 rep("Mesenchymal", 7),
                                                                 rep("Endothelial", 5),rep("Immune", 20)))
split$class = factor(split$class, levels = c("Epithelial","Mesenchymal","Endothelial","Immune"))

#Draw heatmap (Fig. 2a)
Heatmap(startrac.dist, name = "Ro/e", col = colorRamp2(c(min(startrac.dist), 1, max(startrac.dist)), c('#5A8FCA','#F2F2F0','#E31A1C')  ),
        cluster_columns = F, cluster_rows = F, border = T, row_names_side = "left",column_split = split$class )
