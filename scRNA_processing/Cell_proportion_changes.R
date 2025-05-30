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

#Draw heatmap
Heatmap(startrac.dist, name = "Ro/e", col = colorRamp2(c(min(startrac.dist), 1, max(startrac.dist)), c('#5A8FCA','#F2F2F0','#E31A1C')  ),
        cluster_columns = F, cluster_rows = F, border = T, row_names_side = "left",column_split = split$class )

#Calculate cell proportion in each sample
#Epithelial cells
sample_celltype_ratio_epi = lapply(unique(epi_meta$sample), function(sample_id) {
  subset_data = subset(epi_meta, subset = sample == sample_id) 
  celltype_counts = table(subset_data$celltype)
  celltype_ratio_df = data.frame(sample = sample_id, proportion = prop.table(celltype_counts))
  colnames(celltype_ratio_df) = c("sample","celltype","ratio")
  return(celltype_ratio_df)
})
celltype_ratio_epi = do.call(rbind, sample_celltype_ratio_epi)
celltype_ratio_df_epi = dcast(celltype_ratio_epi, sample~celltype, value.var = "ratio")
rownames(celltype_ratio_df_epi) = celltype_ratio_df_epi$sample
celltype_ratio_df_epi$sample = NULL

#Mesenchymal cells
sample_celltype_ratio_mes = lapply(unique(mes_meta$sample), function(sample_id) {
  subset_data = subset(mes_meta, subset = sample == sample_id) 
  celltype_counts = table(subset_data$celltype)
  celltype_ratio_df = data.frame(sample = sample_id, proportion = prop.table(celltype_counts))
  colnames(celltype_ratio_df) = c("sample","celltype","ratio")
  return(celltype_ratio_df)
})
celltype_ratio_mes = do.call(rbind, sample_celltype_ratio_mes)
celltype_ratio_df_mes = dcast(celltype_ratio_mes, sample~celltype, value.var = "ratio")
rownames(celltype_ratio_df_mes) = celltype_ratio_df_mes$sample
celltype_ratio_df_mes$sample = NULL

#Endothelial cells
sample_celltype_ratio_endo = lapply(unique(endo_meta$sample), function(sample_id) {
  subset_data = subset(endo_meta, subset = sample == sample_id) 
  celltype_counts = table(subset_data$celltype)
  celltype_ratio_df = data.frame(sample = sample_id, proportion = prop.table(celltype_counts))
  colnames(celltype_ratio_df) = c("sample","celltype","ratio")
  return(celltype_ratio_df)
})
celltype_ratio_endo = do.call(rbind, sample_celltype_ratio_endo)
celltype_ratio_df_endo = dcast(celltype_ratio_endo, sample~celltype, value.var = "ratio")
rownames(celltype_ratio_df_endo) = celltype_ratio_df_endo$sample
celltype_ratio_df_endo$sample = NULL

#Immune cells
sample_celltype_ratio_immune = lapply(unique(immune_meta$sample), function(sample_id) {
  subset_data = subset(immune_meta, subset = sample == sample_id) 
  celltype_counts = table(subset_data$celltype)
  celltype_ratio_df = data.frame(sample = sample_id, proportion = prop.table(celltype_counts))
  colnames(celltype_ratio_df) = c("sample","celltype","ratio")
  return(celltype_ratio_df)
})
celltype_ratio_immune = do.call(rbind, sample_celltype_ratio_immune)
celltype_ratio_df_immune = dcast(celltype_ratio_immune, sample~celltype, value.var = "ratio")
rownames(celltype_ratio_df_immune) = celltype_ratio_df_immune$sample
celltype_ratio_df_immune$sample = NULL

#Combine cell proportion from each cell class
celltype_ratio_df =  cbind(celltype_ratio_df_epi,celltype_ratio_df_mes,celltype_ratio_df_endo, celltype_ratio_df_immune)
celltype_ratio_df$sample = rownames(celltype_ratio_df)
celltype_ratio_df[is.na(celltype_ratio_df)] = 0

#load clinical information
newmeta = read.csv("tables/meta101.csv", row.names = 1, check.names = FALSE)
newmeta = merge(celltype_ratio_df, newmeta, by = "sample", all = TRUE)

#Calculate correlation between cell proportion and clinical markers
cov_index = list()
index = c("Age","BMI","FEV1%pred","FEV1/FVC(%)","LAA","AWT-PI10")
for (x in index) {
    cov_index[[x]] = data.frame(group = x, celltype = celltype)
    for (cell in celltype)  {
    correlation = cor.test(newmeta[[x]], newmeta[[cell]], method = "pearson")
    p = correlation$p.value
    cor = correlation$estimate
    cov_index[[x]]$p_value[which(cov_index[[x]]$celltype == cell)] = p
    cov_index[[x]]$R[which(cov_index[[x]]$celltype == cell)] = cor
    }
}
cov = do.call(rbind, cov_index)
cov$`-log10(p_value)` = -log10(cov$p_value)

cov$cellclass = "Immune"
cov$cellclass[which(cov$celltype %in% c("Immature AT1","Mature AT1","AT2","TRB Secretory","PreTB Secretory","Basal","Goblet","Ciliated","PNEC"))] = "Epithelial"
cov$cellclass[which(cov$celltype %in% c("Adventitial fibroblast","Alveolar fibroblast","Myofibroblast","Fibromyocyte","Pericyte","ASM","VSM"))] = "Mesenchymal"
cov$cellclass[which(cov$celltype %in% c("Lymphatic","Aerocyte","gCap","Arterial","Venous"))] = "Endothelial"
cov$cellclass = factor(cov$cellclass, levels = c("Epithelial","Mesenchymal","Endothelial","Immune"))
cov$celltype = factor(cov$celltype, levels = c("Immature AT1","Mature AT1","AT2","TRB Secretory","PreTB Secretory","Basal","Goblet","Ciliated","PNEC",
                                               "Adventitial fibroblast","Alveolar fibroblast","Myofibroblast","Fibromyocyte","Pericyte","ASM","VSM",
                                               "Lymphatic","Aerocyte","gCap","Arterial","Venous","B","Plasma","Treg","CD4+ T","CD8+ T","Proliferating T","NKT","NK","cDC1","cDC2","Migratory DC",
                                               "pDC","Classical monocyte","Non-classical monocyte","Alveolar macrophage","Monocyte-derived macrophage","Interstitial macrophage",
                                               "Mast","Neutrophil","Basophil"))
cov$group = factor(cov$group, levels = c("AWT-PI10","LAA","FEV1/FVC(%)","FEV1%pred","BMI","Age"))

ggplot(data = cov, aes(x = celltype, y = group, size = `-log10(p_value)`, col = R)) +
 geom_point() + theme_bw() + xlab("") + ylab("") +  
 facet_grid(. ~ cellclass, scales = "free_x", space = "free_x") +
 scale_color_gradient2(high="#E31A1C",mid = "lightgrey",low ="#5A8FCA", midpoint = 0) +
 theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
       axis.text.y = element_text(color = "black"),       
       panel.grid.major = element_blank(), panel.grid.minor =element_blank())
