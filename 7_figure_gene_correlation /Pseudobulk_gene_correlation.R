library(Seurat)
library(reshape2)
library(ComplexHeatmap)
library(tidyverse)

#load data from Seurat object and clinical information
COPD = readRDS("processed_data_R/COPD_celltype_old_all_counts.rds")
meta = read.csv("./tables/COPD_meta101.csv", row.names = 1, check.names = F)

#calculate pseudobulk gene expression
pseudobulk = AggregateExpression(COPD, assays = "RNA", return.seurat = T, group.by = c("celltype","sample") )
pseudobulk = NormalizeData(pseudobulk, normalization.method = "LogNormalize", scale.factor = 1e6 )

#Add clinical information
pseudobulk@meta.data$celltype = pseudobulk@meta.data$orig.ident
pseudobulk@meta.data$orig.ident = rownames(pseudobulk@meta.data)
pseudobulk@meta.data$sample = sub(".*_", "", rownames(pseudobulk@meta.data))
pseudobulk@meta.data = merge(pseudobulk@meta.data, meta, by = "sample", all = TRUE)
rownames(pseudobulk@meta.data) = pseudobulk@meta.data$orig.ident

#calculate cell numbers in pseudobulk samples
COPD@meta.data$celltype_sample = paste0(COPD@meta.data$celltype,"_",COPD@meta.data$sample)
cellnum = as.data.frame(table(COPD@meta.data$celltype_sample))
colnames(cellnum) = c("group","cellnum")
rownames(cellnum) = cellnum$group
cellnum$group = NULL
pseudobulk = AddMetaData(pseudobulk, metadata = cellnum)

#filter cell numbers less than 10 in each celltype
pseudobulk = subset(pseudobulk, subset = cellnum >= 10)

#calculate gene expression in each cell type
pseudobulk_celltype = AggregateExpression(pseudobulk, assays = "RNA", return.seurat = T, group.by = c("celltype") )
exp_type = GetAssayData(pseudobulk_celltype, assay = "RNA", slot = "data")

#filter out genes expressed in over 20% of cells for each cell type
COPD = readRDS("/datf/mazhuo/jupyter_notebook/COPD/processed_data_R/COPD_celltype.rds")
Idents(COPD) = COPD@meta.data$celltype
genes_select_FEV1 = list()
for (celltype in unique(COPD@meta.data$celltype)) {
  cells_in_type = WhichCells(COPD, idents = celltype)
  percent_expressed = rowSums(COPD@assays$RNA@counts[, cells_in_type] > 0) / length(cells_in_type)
  genes_select_FEV1[[celltype]] = names(percent_expressed[percent_expressed > 0.2])
}

#filter out celltypes less than 50 samples
celltype_num = as.data.frame(table(pseudobulk@meta.data$celltype))
celltype_select = celltype_num$Var1[which(celltype_num$Freq >= 50)]

#calculate Pearson corelation between gene expression and FEV1%pred in each cell type
gene_cor_FEV1 = list()
for (x in celltype_select) {
    cell = subset(pseudobulk, subset = celltype == x)   
    exp = as.data.frame(t(GetAssayData(cell, assay = "RNA", slot = "data")))
    exp$sample = rownames(exp)
    
    cell@meta.data$sample = cell@meta.data$orig.ident
    newmeta = merge(exp, cell@meta.data, by = "sample", all = TRUE)
    rownames(newmeta) = newmeta$sample
    
    gene_cor_FEV1[[x]] = as.data.frame(exp_type[,x])
    colnames(gene_cor_FEV1[[x]]) = "mean_exp"
    gene_cor_FEV1[[x]]$gene = rownames(gene_cor_FEV1[[x]])
    gene_cor_FEV1[[x]] = gene_cor_FEV1[[x]][which(gene_cor_FEV1[[x]]$mean_exp > 1),]
    gene_cor_FEV1[[x]]$celltype = x
    gene_cor_FEV1[[x]]$R = 0
    gene_cor_FEV1[[x]]$p = 0
    
    allgene = intersect(rownames(gene_cor_FEV1[[x]]), genes_select_FEV1[[x]])
    for (gene in allgene) {
    correlation = cor.test(newmeta$`FEV1%pred`, newmeta[[gene]], method = "pearson")
    p = correlation$p.value
    cor = correlation$estimate
    gene_cor_FEV1[[x]]$R[which(gene_cor_FEV1[[x]]$gene == gene)] = cor
    gene_cor_FEV1[[x]]$p[which(gene_cor_FEV1[[x]]$gene == gene)] = p
        }  
    gene_cor_FEV1[[x]] = gene_cor_FEV1[[x]][which(gene_cor_FEV1[[x]]$R != 0 ),]
   }

gene_cor_FEV1_all = do.call(rbind, gene_cor_FEV1)
rownames(gene_cor_FEV1_all) = NULL   
gene_cor_FEV1_all$R_abs = abs(gene_cor_FEV1_all$R)
gene_cor_FEV1_filter = gene_cor_FEV1_all[which(gene_cor_FEV1_all$R_abs > 0.3  ),]

write.csv(gene_cor_FEV1_all, file = "./tables/gene_cor_FEV1_all.csv")
write.csv(gene_cor_FEV1_filter, file = "./tables/gene_cor_FEV1_filter.csv")

#Fig. 5a
gene_cor_FEV1_all$celltype = factor(gene_cor_FEV1_all$celltype, levels = c("Immature AT1","Mature AT1","AT2","TRB Secretory","PreTB Secretory","Basal","Goblet","Ciliated","PNEC",
                 "Adventitial fibroblast","Alveolar fibroblast","Myofibroblast","Fibromyocyte","Pericyte","ASM",
                 "VSM","Mesothelial","Lymphatic","Aerocyte","gCap","Arterial","Venous","B","Plasma","Treg","CD4+ T",
                 "CD8+ T","Proliferating T","NKT","NK","ILC","cDC1","cDC2","Migratory DC","pDC","Classical monocyte","Non-classical monocyte",
                 "Alveolar macrophage","Monocyte-derived macrophage","Interstitial macrophage","Mast","Neutrophil"))
col = setNames(c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896","#9467bd", "#c5b0d5", 
                     "#8c564b", "#c49c94", "#e377c2", "#f786d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d","#17becf","#9edae5", 
                     "#6baed6", "#fd8d3c", "#74c476", "#9e9ac7", "#969696", "#5254a3", "#8ca252", "#bd9e39", "#ad494a", "#a55194",
                     "#FF410D", "#6EE2FF", "#F7C530", "#95cc5e", "#d0dfe6", "#f79d1e", "#748aa6", "#cc0c00", "#5c88da", "#84bd00",
                     "#ffcd00", "#7c878e"),
               c('PNEC','Basal','TRB Secretory','AT2','PreTB Secretory','Goblet','Ciliated','Immature AT1','AT1',
                 'Mesothelial','VSM','ASM','Pericyte','Myofibroblast','Fibromyocyte','Alveolar fibroblast','Adventitial fibroblast',
                 'Lymphatic','gCap','Aerocyte','Arterial','Venous','Migratory DC','cDC1','pDC','ILC','Mast','Plasma','B',
                 'Interstitial macrophage','Monocyte-derived macrophage','Alveolar macrophage','cDC2','Classical monocyte',
                 'Non-classical monocyte','Neutrophil','CD4+ T','Treg','CD8+ T','Proliferating T','NKT','NK'))
label_col = subset(col, names(col) %in% celltype_select)

ggplot(gene_cor_FEV1_all, aes(x = celltype, y = R_abs)) + geom_boxplot(width = 1, outlier.size = 0.1, fill = label_col) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor =element_blank(),axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) ) + scale_color_manual(values = label_col) 

#Fig.5b
celltype_show = "AT2"
gene_show = c("NOP56","RPF2","GPR108","STAT6","TMOD3","LIMA1","CAST","NFKB1")
gene_cor_FEV1_show = gene_cor_FEV1[[celltype_show]]
gene_cor_FEV1_show = gene_cor_FEV1_show %>% arrange(-R)
gene_cor_FEV1_show$rank = 1:nrow(gene_cor_FEV1_show)
gene_cor_FEV1_show$R_abs = abs(gene_cor_FEV1_show$R) 
gene_cor_FEV1_show$group = "|R| <= 0.3"
gene_cor_FEV1_show$group[which(gene_cor_FEV1_show$R > 0.3)] = "R > 0.3"
gene_cor_FEV1_show$group[which(gene_cor_FEV1_show$R < -0.3)] = "R < -0.3"
gene_cor_FEV1_show$label = ""
gene_cor_FEV1_show$label[which(gene_cor_FEV1_show$gene %in% gene_show)] = gene_cor_FEV1_show$gene[which(gene_cor_FEV1_AT2$gene %in% gene_show)]

group_col = setNames(c("grey","#E31A1C","#5A8FCA"),c("|R| <= 0.3","R > 0.3","R < -0.3") )
ggplot(gene_cor_FEV1_show, aes(x = rank, y = R)) + geom_point(size = 1, aes(color = group)) + theme_classic() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) +
   ggtitle(celltype_show) + scale_color_manual(values = group_col) + 
   geom_label_repel(aes(label = label), max.overlaps = 50, min.segment.length = 0, col = "black", label.size = 0.1 )

#Fig. 5c
gene_cor_FEV1_filter_pos = gene_cor_FEV1_filter[which(gene_cor_FEV1_filter$R > 0 ),]
gene_cor_FEV1_filter_pos_num = as.data.frame(table(gene_cor_FEV1_filter_pos$celltype))
colnames(gene_cor_FEV1_filter_pos_num) = c("celltype","genenum")
gene_cor_FEV1_filter_pos_num$group = "Positively"

gene_cor_FEV1_filter_neg = gene_cor_FEV1_filter[which(gene_cor_FEV1_show$R < 0 ),]
gene_cor_FEV1_filter_neg_num = as.data.frame(table(gene_cor_FEV1_filter_neg$celltype))
colnames(gene_cor_FEV1_filter_neg_num) = c("celltype","genenum")
gene_cor_FEV1_filter_neg_num$group = "Negatively"

gene_cor_FEV1_show_all_num = rbind(gene_cor_FEV1_show_pos_num, gene_cor_FEV1_show_neg_num)
group_col = setNames(c('#E31A1C','#5A8FCA'), c("Positively","Negatively") )

ggplot(gene_cor_FEV1_filter_all_num, aes(x = celltype, y = genenum, fill = group)) +
  geom_bar(stat = "identity", position = "stack") + scale_fill_manual(values = group_col) +
  labs(title = "Genes associated with FEV1%predicted", x = "", y = "# gene") +  theme_classic() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor =element_blank(), plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, h = 1)) + 
  scale_x_discrete(limits = celltype_select)

#Fig .5e-f
celltype_show = "AT2"
gene_show = "NFKB1"
cell = subset(pseudobulk, subset = celltype == celltype_show)
matrix = cell@assays$RNA@data
colnames(matrix) = gsub("_.*", "", colnames(matrix) )
matrix = as.data.frame(t(matrix))
matrix$sample = rownames(matrix)
matrix = left_join(matrix, meta, by = "sample")

group_col = setNames(c("#51a246","#3194cc","#794191","#684b3b"),
                     c("HC","HS","Mild COPD","Severe COPD") )
ggplot(matrix, aes(x = `FEV1%pred`, y = !!sym(gene_show)  )  ) + 
         geom_point(aes(color = new_group)) + geom_smooth(method = "lm") + scale_color_manual(values = group_col) +
         theme_classic() +theme(panel.grid.major = element_blank(), panel.grid.minor =element_blank(),
                                plot.title = element_text(hjust = 0.5)) + 
         stat_cor(method = "pearson",size = 5) +
         ylab("Gene expression") + xlab("FEV1%pred") + ggtitle(paste0(gene_show, "in", celltype_show) )
