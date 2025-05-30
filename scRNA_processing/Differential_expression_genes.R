library(Seurat)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(ComplexHeatmap)
library(colorRamp2)

#Load Seurat object
COPD = readRDS("processed_data_R/COPD_celltype_counts_0808_v3.rds")
COPD = NormalizeData(COPD)
Idents(COPD) = COPD@meta.data$Group
HCHS = subset(COPD, subset = Group %in% c("HC","HS") )
HSG12 = subset(COPD, subset = Group %in% c("HS","Mild COPD") )
G1234 = subset(COPD, subset = Group %in% c("Mild COPD","Severe COPD") )
rm(COPD)

celltype = c("Immature AT1","Mature AT1","AT2","TRB Secretory","PreTB Secretory","Basal","Goblet","Ciliated","PNEC",
                 "Adventitial fibroblast","Alveolar fibroblast","Myofibroblast","Fibromyocyte","Pericyte","ASM",
                 "VSM","Mesothelial","Lymphatic","Aerocyte","gCap","Arterial","Venous","B","Plasma","Treg","CD4+ T",
                 "CD8+ T","Proliferating T","NKT","NK","Basophil","cDC1","cDC2","Migratory DC","pDC","Classical monocyte","Non-classical monocyte",
                 "Alveolar macrophage","Monocyte-derived macrophage","Interstitial macrophage","Mast","Neutrophil")

#DE analysis between HC and HS group
DEG_wilcox_HCHS = list()
celltype = unique(as.character(HCHS@meta.data$celltype))
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
DEG_wilcox_HCHS = DEG_wilcox_HCHS[which(DEG_wilcox_HCHS$p_val_adj < 0.05),]
DEG_wilcox_HCHS = DEG_wilcox_HCHS[which(DEG_wilcox_HCHS$pct.1 > 0.2),]
write.csv(DEG_wilcox_HCHS, file = "tables/COPD_DEG_wilcox_HCHS.csv" )

#DE analysis between HS and GOLD1-2 group
DEG_wilcox_HSG12 = list()
celltype = unique(as.character(HSG12@meta.data$celltype))
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
DEG_wilcox_HSG12 = DEG_wilcox_HSG12[which(DEG_wilcox_HSG12$p_val_adj < 0.05),]
DEG_wilcox_HSG12 = DEG_wilcox_HSG12[which(DEG_wilcox_HSG12$pct.1 > 0.2),]
write.csv(DEG_wilcox_HSG12, file = "tables/COPD_DEG_wilcox_HSG12.csv" )

#DE analysis between GOLD1-2 and GOLD3-4 group
DEG_wilcox_G1234 = list()
celltype = unique(as.character(G1234@meta.data$celltype))
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
DEG_wilcox_G1234 = DEG_wilcox_G1234[which(DEG_wilcox_G1234$p_val_adj < 0.05),]
DEG_wilcox_G1234 = DEG_wilcox_G1234[which(DEG_wilcox_G1234$pct.1 > 0.2),]
write.csv(DEG_wilcox_G1234, file = "tables/COPD_DEG_wilcox_G1234.csv" )

#Calculate up regulated DEG numbers
DEG_HCHS_num = DEG_wilcox_HCHS[which(DEG_wilcox_HCHS$avg_log2FC < 0 ),]
DEG_HCHS_num = as.data.frame(table(DEG_HCHS_num$celltype))
colnames(DEG_HCHS_num) = c("celltype","genenum")
DEG_HCHS_num$group = "Pre-change"

DEG_HSG12_num = DEG_wilcox_HSG12[which(DEG_wilcox_HSG12$avg_log2FC < 0 ),]
DEG_HSG12_num = as.data.frame(table(DEG_HSG12_num$celltype))
colnames(DEG_HSG12_num) = c("celltype","genenum")
DEG_HSG12_num$group = "Early change"

DEG_G1234_num = DEG_wilcox_G1234[which(DEG_wilcox_G1234$avg_log2FC < 0 ),]
DEG_G1234_num = as.data.frame(table(DEG_G1234_num$celltype))
colnames(DEG_G1234_num) = c("celltype","genenum")
DEG_G1234_num$group = "Late change"

DEG_group_up = rbind(DEG_HCHS_num, DEG_HSG12_num, DEG_G1234_num)
DEG_group_up$group = factor(DEG_group_up$group, levels = c("Pre-change","Early change","Late change") )

DEG_up_num = aggregate(genenum ~ celltype, data = DEG_group_up, sum)
DEG_up_num = setNames(DEG_up_num$genenum,DEG_up_num$celltype)

DEG_group_up_num = aggregate(genenum ~ group, data = DEG_group_up, sum)
DEG_group_up_num = setNames(DEG_group_up_num$genenum,DEG_group_up_num$group)

#Calculate down regulated DEG numbers
DEG_HCHS_num = DEG_wilcox_HCHS[which(DEG_wilcox_HCHS$avg_log2FC > 0 ),]
DEG_HCHS_num = as.data.frame(table(DEG_HCHS_num$celltype))
colnames(DEG_HCHS_num) = c("celltype","genenum")
DEG_HCHS_num$group = "Pre-change"

DEG_HSG12_num = DEG_wilcox_HSG12[which(DEG_wilcox_HSG12$avg_log2FC > 0 ),]
DEG_HSG12_num = as.data.frame(table(DEG_HSG12_num$celltype))
colnames(DEG_HSG12_num) = c("celltype","genenum")
DEG_HSG12_num$group = "Early change"

DEG_G1234_num = DEG_wilcox_G1234[which(DEG_wilcox_G1234$avg_log2FC > 0 ),]
DEG_G1234_num = as.data.frame(table(DEG_G1234_num$celltype))
colnames(DEG_G1234_num) = c("celltype","genenum")
DEG_G1234_num$group = "Late change"

DEG_group_down = rbind(DEG_HCHS_num, DEG_HSG12_num, DEG_G1234_num)
DEG_group_down$group = factor(DEG_group_down$group, levels = c("Pre-change","Early change","Late change"))

DEG_down_num = aggregate(genenum ~ celltype, data = DEG_group_down, sum)
DEG_down_num = setNames(DEG_down_num$genenum,DEG_down_num$celltype)

DEG_group_down_num = aggregate(genenum ~ group, data = DEG_group_down, sum)
DEG_group_down_num = setNames(DEG_group_down_num$genenum,DEG_group_down_num$group)

#Draw Heatmap of DEG numbers
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

col = setNames(c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896","#9467bd", "#c5b0d5", 
                     "#8c564b", "#c49c94", "#e377c2", "#f786d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d","#17becf","#9edae5", 
                     "#6baed6", "#fd8d3c", "#74c476", "#9e9ac7", "#969696", "#5254a3", "#8ca252", "#bd9e39", "#ad494a", "#a55194",
                     "#FF410D", "#6EE2FF", "#F7C530", "#95cc5e", "#d0dfe6", "#f79d1e", "#748aa6", "#cc0c00", "#5c88da", "#84bd00",
                     "#ffcd00", "#7c878e"),
               c('PNEC','Basal','TRB Secretory','AT2','PreTB Secretory','Goblet','Ciliated','Immature AT1','Mature AT1',
                 'Mesothelial','VSM','ASM','Pericyte','Myofibroblast','Fibromyocyte','Alveolar fibroblast','Adventitial fibroblast',
                 'Lymphatic','gCap','Aerocyte','Arterial','Venous','Migratory DC','cDC1','pDC','Basophil','Mast','Plasma','B',
                 'Interstitial macrophage','Monocyte-derived macrophage','Alveolar macrophage','cDC2','Classical monocyte',
                 'Non-classical monocyte','Neutrophil','CD4+ T','Treg','CD8+ T','Proliferating T','NKT','NK'))

up_col = setNames(rep('#E31A1C',3), c("Pre-change","Early change","Late change"))
down_col = setNames(rep('#5A8FCA',3), c("Pre-change","Early change","Late change"))

split = data.frame(celltype = colnames(DEG_up_matrix), class = c(rep("Epithelial", 9),
                                                                 rep("Mesenchymal", 8),
                                                                 rep("Endothelial", 5),rep("Immune", 20)))
split$class = factor(split$class, levels = c("Epithelial","Mesenchymal","Endothelial","Immune"))

celltype_df = as.data.frame(celltype)
rownames(celltype_df) = celltype_df$celltype
col_anno = HeatmapAnnotation(celltype = celltype_df$celltype, col=list(celltype = col))

ha1 = HeatmapAnnotation(genenum = anno_barplot(DEG_up_num, height = unit(2, "cm"), gp = gpar(fill = col[colnames(DEG_up_matrix)], col = col[colnames(DEG_up_matrix)])) )                     
haa = rowAnnotation(genenum = anno_barplot(DEG_group_up_num, height = unit(2, "cm"), gp = gpar(fill = up_col, col = up_col)))
ht1 = Heatmap(DEG_up_matrix, col = colorRamp2(c(0,800), c('#F2F2F0','#E31A1C')), 
              name = "Upregualted DEGs", row_names_side = "left", 
              cluster_rows = F, cluster_columns = F, border = T,
              top_annotation = ha1, right_annotation = haa,
              bottom_annotation = col_anno,
              heatmap_legend_param = list(legend_direction = "horizontal"),
              column_split = split$class)

ha2 = HeatmapAnnotation(genenum = anno_barplot(DEG_down_num, height = unit(2, "cm"), gp = gpar(fill = col[colnames(DEG_up_matrix)], col = col[colnames(DEG_up_matrix)])) )
hab = rowAnnotation(genenum = anno_barplot(DEG_group_down_num, height = unit(2, "cm"), gp = gpar(fill = down_col, col = down_col) ))
ht2 = Heatmap(DEG_down_matrix, col = colorRamp2(c(0,800), c('#F2F2F0','#5A8FCA')), name = "Downregualted DEGs",
             row_names_side = "left", cluster_rows = F, cluster_columns = F, border = T,
             top_annotation = ha2, right_annotation = hab,
              bottom_annotation = col_anno,
              heatmap_legend_param = list(legend_direction = "horizontal"),
              column_split = split$class )

ht_list = ht1 %v% ht2
draw(ht_list)
