library(Seurat)
library(reshape2)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(colorRamp2)
set.seed(123)

COPD = readRDS("processed_data_R/COPD_celltype_counts_0808_v3.rds")
COPD = NormalizeData(COPD)
epi = subset(COPD, subset = cellclass == "Epithelial")

HCHS = subset(epi, subset = Group %in% c("HC","HS") )
HSG12 = subset(epi, subset = Group %in% c("HS","Mild COPD") )
G1234 = subset(epi, subset = Group %in% c("Mild COPD","Severe COPD") )

#Pre-change
DEG_wilcox_HCHS = list()
celltype = c("Immature AT1","Mature AT1","AT2","TRB Secretory","PreTB Secretory","Basal","Goblet","Ciliated","PNEC")
for (x in celltype) {
    cell = subset(HCHS, subset = celltype == x)
    DEG_wilcox_HCHS[[x]] = FindMarkers(cell, ident.1 = "HC", ident.2 = "HS", logfc.threshold = 0.25)
    DEG_wilcox_HCHS[[x]]$group = "HS"
    DEG_wilcox_HCHS[[x]]$group[which(DEG_wilcox_HCHS[[x]]$avg_log2FC > 0)] = "HC"
    DEG_wilcox_HCHS[[x]]$gene = rownames(DEG_wilcox_HCHS[[x]])
    DEG_wilcox_HCHS[[x]]$celltype = x
}
DEG_wilcox_HCHS = do.call(rbind, DEG_wilcox_HCHS)
rownames(DEG_wilcox_HCHS) = NULL
DEG_wilcox_HCHS$abs_log2FC = abs(DEG_wilcox_HCHS$avg_log2FC)
DEG_wilcox_HCHS = DEG_wilcox_HCHS[which(DEG_wilcox_HCHS$p_val_adj < 0.05 &
                                        DEG_wilcox_HCHS$abs_log2FC > 0.3 &
                                        DEG_wilcox_HCHS$pct.1 > 0.2) ,]
write.csv(DEG_wilcox_HCHS, file = "./tables/COPD_DEG_wilcox_HCHS_epi_filter.csv" )

#Early change
DEG_wilcox_HSG12 = list()
for (x in celltype) {
    cell = subset(HSG12, subset = celltype == x)
    DEG_wilcox_HSG12[[x]] = FindMarkers(cell, ident.1 = "HS", ident.2 = "Mild COPD", logfc.threshold = 0.25)
    DEG_wilcox_HSG12[[x]]$group = "Mild COPD"
    DEG_wilcox_HSG12[[x]]$group[which(DEG_wilcox_HSG12[[x]]$avg_log2FC > 0)] = "HS"
    DEG_wilcox_HSG12[[x]]$gene = rownames(DEG_wilcox_HSG12[[x]])
    DEG_wilcox_HSG12[[x]]$celltype = x
}
DEG_wilcox_HSG12 = do.call(rbind, DEG_wilcox_HSG12)
rownames(DEG_wilcox_HSG12) = NULL
write.csv(DEG_wilcox_HSG12, file = "tables/COPD_DEG_wilcox_HSG12_epi.csv" )
DEG_wilcox_HSG12$abs_log2FC = abs(DEG_wilcox_HSG12$avg_log2FC)
DEG_wilcox_HSG12 = DEG_wilcox_HSG12[which(DEG_wilcox_HSG12$p_val_adj < 0.05 &
                                          DEG_wilcox_HSG12$abs_log2FC > 0.3 &
                                          DEG_wilcox_HSG12$pct.1 > 0.2) ,]
write.csv(DEG_wilcox_HSG12, file = "./tables/COPD_DEG_wilcox_HSG12_epi_filter.csv" )

#Late change
DEG_wilcox_G1234 = list()
for (x in celltype) {
    cell = subset(G1234, subset = celltype == x)
    DEG_wilcox_G1234[[x]] = FindMarkers(cell, ident.1 = "Mild COPD", ident.2 = "Severe COPD", logfc.threshold = 0.25)
    DEG_wilcox_G1234[[x]]$group = "Severe COPD"
    DEG_wilcox_G1234[[x]]$group[which(DEG_wilcox_G1234[[x]]$avg_log2FC > 0)] = "Mild COPD"
    DEG_wilcox_G1234[[x]]$gene = rownames(DEG_wilcox_G1234[[x]])
    DEG_wilcox_G1234[[x]]$celltype = x
}
DEG_wilcox_G1234 = do.call(rbind, DEG_wilcox_G1234)
rownames(DEG_wilcox_G1234) = NULL
write.csv(DEG_wilcox_G1234, file = "tables/COPD_DEG_wilcox_G1234_epi.csv" )
DEG_wilcox_G1234$abs_log2FC = abs(DEG_wilcox_G1234$avg_log2FC)
DEG_wilcox_G1234 = DEG_wilcox_G1234[which(DEG_wilcox_G1234$p_val_adj < 0.05 &
                                         DEG_wilcox_G1234$abs_log2FC > 0.3 &
                                         DEG_wilcox_G1234$pct.1 > 0.2) ,]
write.csv(DEG_wilcox_G1234, file = "./tables/COPD_DEG_wilcox_G1234_epi_filter.csv" )


up_HCHS_num = DEG_wilcox_HCHS[which(DEG_wilcox_HCHS$avg_log2FC < 0 ),]
up_HCHS_num = as.data.frame(table(up_HCHS_num$celltype))
colnames(up_HCHS_num) = c("celltype","genenum")
up_HCHS_num$group = "Pre-change"

up_HSG12_num = DEG_wilcox_HSG12[which(DEG_wilcox_HSG12$avg_log2FC < 0 ),]
up_HSG12_num = as.data.frame(table(up_HSG12_num$celltype))
colnames(up_HSG12_num) = c("celltype","genenum")
up_HSG12_num$group = "Early change"

up_G1234_num = DEG_wilcox_G1234[which(DEG_wilcox_G1234$avg_log2FC < 0 ),]
up_G1234_num = as.data.frame(table(up_G1234_num$celltype))
colnames(up_G1234_num) = c("celltype","genenum")
up_G1234_num$group = "Late change"

DEG_group_up = rbind(up_HCHS_num, up_HSG12_num, up_G1234_num)
DEG_group_up$group = factor(DEG_group_up$group, levels = c("Pre-change","Early change","Late change") )
DEG_up_num = aggregate(genenum ~ celltype, data = DEG_group_up, sum)
DEG_up_num = setNames(DEG_up_num$genenum,DEG_up_num$celltype)
DEG_group_up_num = aggregate(genenum ~ group, data = DEG_group_up, sum)
DEG_group_up_num = setNames(DEG_group_up_num$genenum,DEG_group_up_num$group)

down_HCHS_num = DEG_wilcox_HCHS[which(DEG_wilcox_HCHS$avg_log2FC > 0 ),]
down_HCHS_num = as.data.frame(table(down_HCHS_num$celltype))
colnames(down_HCHS_num) = c("celltype","genenum")
down_HCHS_num$group = "Pre-change"

down_HSG12_num = DEG_wilcox_HSG12[which(DEG_wilcox_HSG12$avg_log2FC > 0 ),]
down_HSG12_num = as.data.frame(table(down_HSG12_num$celltype))
colnames(down_HSG12_num) = c("celltype","genenum")
down_HSG12_num$group = "Early change"

down_G1234_num = DEG_wilcox_G1234[which(DEG_wilcox_G1234$avg_log2FC > 0 ),]
down_G1234_num = as.data.frame(table(down_G1234_num$celltype))
colnames(down_G1234_num) = c("celltype","genenum")
down_G1234_num$group = "Late change"

DEG_group_down = rbind(down_HCHS_num, down_HSG12_num, down_G1234_num)
DEG_group_down$group = factor(DEG_group_down$group, levels = c("Pre-change","Early change","Late change"))
DEG_down_num = aggregate(genenum ~ celltype, data = DEG_group_down, sum)
DEG_down_num = setNames(DEG_down_num$genenum,DEG_down_num$celltype)
DEG_group_down_num = aggregate(genenum ~ group, data = DEG_group_down, sum)
DEG_group_down_num = setNames(DEG_group_down_num$genenum,DEG_group_down_num$group)

DEG_up_matrix = dcast(DEG_group_up, group ~ celltype, value.var = "genenum", drop = F )
DEG_down_matrix = dcast(DEG_group_down, group ~ celltype, value.var = "genenum", drop = F )
DEG_up_matrix[is.na(DEG_up_matrix)] = 0
DEG_down_matrix[is.na(DEG_down_matrix)] = 0
rownames(DEG_up_matrix) = DEG_up_matrix$group
DEG_up_matrix$group = NULL
DEG_up_matrix = as.matrix(DEG_up_matrix)
rownames(DEG_down_matrix) = DEG_down_matrix$group
DEG_down_matrix$group = NULL
DEG_down_matrix = as.matrix(DEG_down_matrix)

epi_col = setNames(c("#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#b3bd22", "#17becf", "#aec7e8"),
                   c("Immature AT1","AT1","AT2","TRB Secretory","PreTB Secretory","Basal","Goblet","Ciliated","PNEC"))
up_col = setNames(rep('#E31A1C',3), c("Pre-change","Early change","Late change"))
down_col = setNames(rep('#5A8FCA',3), c("Pre-change","Early change","Late change"))

celltype_df = as.data.frame(celltype)
rownames(celltype_df) = celltype_df$celltype
col_anno = HeatmapAnnotation(celltype = celltype_df$celltype, col=list(celltype = epi_col))

ha1 = HeatmapAnnotation(genenum = anno_barplot(DEG_up_num, height = unit(2, "cm"), 
                                               gp = gpar(fill = epi_col[colnames(DEG_up_matrix)], 
                                                         col = epi_col[colnames(DEG_up_matrix)])) )                     
haa = rowAnnotation(genenum = anno_barplot(DEG_group_up_num, height = unit(2, "cm"), 
                                           gp = gpar(fill = up_col, col = up_col)))
ht1 = Heatmap(DEG_up_matrix, col = colorRamp2(c(0,300), c('#F2F2F0','#E31A1C')), 
              name = "Upregualted DEGs", row_names_side = "left", 
              cluster_rows = F, cluster_columns = F, border = T,
              top_annotation = ha1, right_annotation = haa,
              bottom_annotation = col_anno,
              heatmap_legend_param = list(legend_direction = "horizontal"))

ha2 = HeatmapAnnotation(genenum = anno_barplot(DEG_down_num, height = unit(2, "cm"), 
                                               gp = gpar(fill = epi_col[colnames(DEG_up_matrix)], 
                                                         col = epi_col[colnames(DEG_up_matrix)])) )
hab = rowAnnotation(genenum = anno_barplot(DEG_group_down_num, height = unit(2, "cm"), 
                                           gp = gpar(fill = down_col, col = down_col) ))
ht2 = Heatmap(DEG_down_matrix, col = colorRamp2(c(0,300), c('#F2F2F0','#5A8FCA')), name = "Downregualted DEGs",
             row_names_side = "left", cluster_rows = F, cluster_columns = F, border = T,
             top_annotation = ha2, right_annotation = hab,
              bottom_annotation = col_anno,
              heatmap_legend_param = list(legend_direction = "horizontal") )

ht_list = ht1 %v% ht2
draw(ht_list)
