library(Seurat)
library(reshape2)
library(ComplexHeatmap)
library(tidyverse)

#load data from Seurat object and clinical information
COPD = readRDS("processed_data_R/COPD_celltype_old_all_counts.rds")
meta = read.csv("./tables/meta101.csv", row.names = 1, check.names = F)

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
COPD = readRDS("/datf/mazhuo/jupyter_notebook/COPD/processed_data_R/COPD_celltype_counts_0808_v3.rds")
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

write.csv(gene_cor_FEV1_all, file = "./tables/gene_cor_FEV1_all_241219.csv")
write.csv(gene_cor_FEV1_filter, file = "./tables/gene_cor_FEV1_filter_241219.csv")
