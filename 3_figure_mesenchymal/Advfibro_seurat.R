library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(EnhancedVolcano)
library(patchwork)
library(reshape2)
library(ggpubr)
library(ComplexHeatmap)
library(colorRamp2)
library(harmony)
library(ggpubr)
library(rstatix)

COPD = readRDS("./processed_data_R/COPD_celltype_umap.rds")
advfibro = subset(COPD, subset = celltype == "Adventitial fibroblast")
advfibro_meta = read.csv("./tables/advfibro_leiden_meta.csv", row.names = 1, check.names = FALSE)
advfibro@meta.data = advfibro_meta
UMAP = as.matrix(advfibro@meta.data[,c("UMAP1","UMAP2")])
colnames(UMAP) = c("umap_1","umap_2")
advfibro@reductions$umap@cell.embeddings = UMAP

advfibro@meta.data$cellsubtype = factor(advfibro@meta.data$cellsubtype,
                                       levels = c("Advfibro_1","Advfibro_2","Advfibro_3") )
Idents(advfibro) = advfibro@meta.data$cellsubtype
adv_col = setNames(c('#2ca02c', '#d62728', '#9467bd'),c("Advfibro_1","Advfibro_2","Advfibro_3"))

DimPlot(advfibro, reduction = "umap", group.by = "cellsubtype", 
        pt.size = 1, label = TRUE) + scale_color_manual(values = adv_col) 

Adv_gene = list("Advfibro_1" = c("RORA","PRKG1","ZNF385D"),
                "Advfibro_2" = c("PTX3","MT1A","TNFAIP6","CCL2","CXCL2","IL6","RUNX1"),
                "Advfibro_3" = c("PI16","MFAP5","SLPI"))
DotPlot(advfibro, features = Adv_gene, group.by = "cellsubtype") + RotatedAxis() + 
     scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0) 

advfibro_marker = FindAllMarkers(advfibro, logfc.threshold = 0.4, only.pos = TRUE)
advfibro_marker = advfibro_marker[which(advfibro_marker$p_val_adj < 0.05),]

advfibro@meta.data$cellsubtype = factor(advfibro@meta.data$cellsubtype, 
                                     levels = c("Advfibro_1","Advfibro_2","Advfibro_3") )
advfibro@meta.data = advfibro@meta.data[order(advfibro@meta.data$cellsubtype), ]
advfibro_meta = subset(advfibro@meta.data, select = "cellsubtype")

matrix = as.matrix(advfibro@assays$RNA@data[advfibro_marker$gene,rownames(advfibro_meta)])
matrix = t(apply(matrix,1,function(x){(x-mean(x))/sd(x)}))
col_anno = HeatmapAnnotation(group = advfibro_meta$cellsubtype,
                             col = list(group = adv_col))

Heatmap(matrix, col = colorRamp2(seq(from=-2,to=2,length=11),colorRampPalette(c("#5A8FCA", "white", "#E31A1C"))(11)),
        name = "Expression", cluster_columns = FALSE, cluster_rows = FALSE, use_raster = TRUE, 
        top_annotation = col_anno,
        column_split = c(rep("Advfibro_1", 16991), rep("Advfibro_2", 13596), rep("Advfibro_3", 12211)),
        row_split = c(rep("Advfibro_1", 124), rep("Advfibro_2", 349), rep("Advfibro_3", 197)),     
       show_row_names = FALSE,  show_column_names = FALSE, border = TRUE ) 

Gene_ID =bitr(advfibro_marker$gene, fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")
advfibro_marker = merge(Gene_ID,advfibro_marker,by.x='SYMBOL',by.y='gene')
advfibro_marker_GO = compareCluster(
  ENTREZID ~ cluster, 
  data=advfibro_marker , 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
dotplot(advfibro_marker_GO, showCategory = 4) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1))

sample_celltype_ratio = lapply(unique(advfibro_meta$sample), function(sample_id) {
  subset_data = subset(advfibro_meta, subset = sample == sample_id) 
  celltype_counts = table(subset_data$cellsubtype)
  celltype_ratio_df = data.frame(sample = sample_id, proportion = prop.table(celltype_counts))
  colnames(celltype_ratio_df) = c("sample","cellsubtype","ratio")
  return(celltype_ratio_df)
})
celltype_ratio = do.call(rbind, sample_celltype_ratio)
celltype_ratio_df = dcast(celltype_ratio, sample~cellsubtype, value.var = "ratio")
rownames(celltype_ratio_df) = celltype_ratio_df$sample
celltype_ratio_df[is.na(celltype_ratio_df)] = 0

cellnum = as.data.frame(table(advfibro_meta$sample))
sample = cellnum[which(cellnum$Freq >= 50),]$Var1
newmeta = read.csv("tables/COPD_meta101.csv", row.names = 1, check.names = FALSE)
newmeta$sample = rownames(newmeta)
celltype_ratio_df = left_join(celltype_ratio_df, newmeta, by = "sample")
celltype_ratio_df = celltype_ratio_df[which(celltype_ratio_df$sample %in% sample),]

group_col = setNames(c("#51a246","#3194cc","#794191","#684b3b"),
                     c("HC","HS","Mild COPD","Severe COPD") )

plot = list()
for (x in celltype) {
    p = ggplot(celltype_ratio_df, aes(x = `FEV1%pred`, y = !!sym(x))  ) + 
           geom_point(aes(color = Group)) + geom_smooth(method = "lm") + scale_color_manual(values = group_col) +
           theme_classic() +theme(panel.grid.major = element_blank(), panel.grid.minor =element_blank(),
                                  plot.title = element_text(hjust = 0.5),
                                 legend.position = "none") + 
           stat_cor(method = "pearson",size = 5) +
           ylab("Proportion") + xlab("FEV1%pred") + ggtitle(x)
    plot[[x]] = p
    }
wrap_plots(plot, ncol = 3)
