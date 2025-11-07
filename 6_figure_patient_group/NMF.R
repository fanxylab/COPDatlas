library(NMF)
library(ComplexHeatmap)
library(reshape2)
library(tidyverse)
library(dplyr)
library(colorRamp2)
library(ggpubr)
library(patchwork)
library(rstatix)
library(ggplot2)
library(Seurat)
library(epitools)

#load cell ratio information
cellratio_df = read.csv('./tables/COPD_cellratio_new.csv', 
                        row.names = 1, check.names = FALSE)
celltype = c('B', 'Plasma', 'CD4+ T', 'Treg', 'CD8+ T', 'NKT', 'NK', 
             'cDC1', 'cDC2', 'Migratory DC', 'pDC', 'Classical monocyte', 
             'Non-classical monocyte', 'Alveolar macrophage', 'Monocyte-derived macrophage', 
             'Interstitial macrophage', 'Mast', 'Neutrophil')

meta = read.csv("./tables/COPD_meta101.csv", row.names = 1, check.names = FALSE)
meta$sample = rownames(meta)
sample = meta$sample[which(meta$Disease == "COPD")]
meta = meta[sample,]

c_df = cellratio_df[sample, celltype]
c_df = t(as.data.frame(c_df))
ranks = 2:10
estim.coad = nmf(c_df, ranks, nrun=20, seed = 32)
plot(estim.coad)

nmf.rank = nmf(c_df, rank = 2, nrun= 20,  seed = 32)

NMF_group = as.data.frame(predict(nmf.rank))
colnames(NMF_group) = "NMF_group"
NMF_group$NMF_group = paste0("T",as.character(NMF_group$NMF_group))
NMF_group$sample = rownames(NMF_group)
NMF_group = left_join(meta, NMF_group, by = "sample")
rownames(NMF_group) = NMF_group$sample

basis_celltype = basis(nmf.rank)
colnames(basis_celltype) = c("F1","F2")
basis_celltype_scale = apply(basis_celltype, 2, function(x) {
  (x - min(x)) / (max(x) - min(x))
})

Heatmap(t(basis_celltype_scale),  name = "Loading", 
              cluster_rows = FALSE, border = TRUE, row_names_side = c("left"),
              col = colorRamp2(c(0, 0.5, 1), c('#5A8FCA','#F2F2F0','#E31A1C')), 
              column_title = "Cell type")

coef_patients = coef(nmf.rank)
rownames(coef_patients) = c("F1","F2")
coef_patients_scale = apply(coef_patients, 1, function(x) {
  (x - min(x)) / (max(x) - min(x))
})

ha = HeatmapAnnotation(NMF_group = NMF_group$NMF_group, 
                       col = list(NMF_group = c("T1" = "#1f77b4","T2" = "#ff7f0e")) )
Heatmap(t(coef_patients_scale),  name = "Loading", cluster_rows = FALSE, border = TRUE,
              col = colorRamp2(c(0, 0.5, 1), c('#5A8FCA','#F2F2F0','#E31A1C')),
              row_names_side = c("left"),  
              #show_column_names = FALSE,
              top_annotation = ha, column_title = "Patient")

index = extractFeatures(nmf.rank,"max") 
sig.order = unlist(index)
NMF.Exp.rank = c_df[sig.order,]
NMF.Exp.rank = na.omit(NMF.Exp.rank) 

row.group = c()
for(each in rownames(NMF.Exp.rank)){
  if(each %in% rownames(c_df)[index[[1]]]){
    row.group = c(row.group, "F1")
  }else if(each %in% rownames(c_df)[index[[2]]]){
    row.group = c(row.group, "F2")
  }
}

col.group = predict(nmf.rank)
col.group = factor(col.group, levels = c("1","2"))
ratio_z = t(scale(cellratio_df))
ratio_z = ratio_z[sig.order,]

hb = HeatmapAnnotation(Stage = factor(NMF_group$Group, levels = c("Mild COPD","Severe COPD")),
                       AECOPD_group = factor(NMF_group$AECOPD_group, 
                                             levels = c("Not frequent AE","Frequent AE")),
                       AECOPD = factor(NMF_group$AECOPD, levels = c("0","1","2","3")),
                       NMF_group = factor(NMF_group$NMF_group, levels = c("T1", "T2") ),
                       col = list(Stage = c("Mild COPD" = "#7967BD","Severe COPD" = "#8c564b"),
                                  AECOPD_group = c("Not frequent AE" = "#5A8FCA",
                                                   "Frequent AE" = "#E31A1C"),
                                  AECOPD = c("0" = "#C6DBEF","1" = "#9ECAE1",
                                             "2" = "#6BAED6","3" = "#3182DD"),
                                  NMF_group = c("T1" = "#1f77b4","T2" = "#ff7f0e") ))

#Fig. 4e
Heatmap(ratio_z, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
        col = colorRamp2(c(-2, 0, 2), c('#5A8FCA','#F2F2F0','#E31A1C')  ), border = TRUE,
        row_split = row.group, name = "Relative cell proportion in immune cells",
        column_split = col.group, top_annotation = hb)


