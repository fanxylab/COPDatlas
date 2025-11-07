library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(colorRamp2)

#load cell metadata
epi_meta = read.csv("./tables/epi_meta_umap.csv")
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

#load clinical information
newmeta = read.csv("tables/COPD_meta101.csv", row.names = 1, check.names = FALSE)
newmeta = merge(celltype_ratio_df, newmeta, by = "sample", all = TRUE)

#Calculate correlation between cell proportion and clinical markers
cov_index = list()
index = c("Age","BMI","FEV1%pred","FEV1/FVC(%)","LAA","AWT-PI10","Pack_years","AECOPD",
          "CAT_score","mMRC_score")
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

cov_ht = dcast(cov, celltype~group, value.var = "R", drop = F)
rownames(cov_ht) = cov_ht$celltype
cov_ht$celltype = NULL

p_ht = dcast(cov, celltype~group, value.var = "p_value", drop = F)
rownames(p_ht) = p_ht$celltype
p_ht$celltype = NULL

significant_label <- ifelse(p_ht < 0.01, "**",
                           ifelse(p_ht < 0.05, "*",
                                 ifelse(p_ht < 0.1, "#", "")))

Heatmap(cov_ht, col = colorRamp2(c(-0.6,0,0.6), c('#5A8FCA','#F2F2F0','#E31A1C') ), 
             name = "Correlation", cluster_rows = F, cluster_columns = F, 
             border = T, row_names_side = c("left"),
             cell_fun = function(j, i, x, y, width, height, fill) {
             grid.text(significant_label[i,j], x, y, gp = gpar(fontsize = 10))
