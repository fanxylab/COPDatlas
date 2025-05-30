
#load data from Seurat object
COPD = readRDS("processed_data_R/COPD_celltype_old_all_counts.rds")

#calculate pseudobulk gene expression
pseudobulk = AggregateExpression(COPD, assays = "RNA", return.seurat = T, group.by = c("celltype","sample") )
pseudobulk = NormalizeData(pseudobulk, normalization.method = "LogNormalize", scale.factor = 1e6 )

#calculate cell numbers in pseudobulk samples
COPD@meta.data$celltype_sample = paste0(COPD@meta.data$celltype,"_",COPD@meta.data$sample)
cellnum = as.data.frame(table(COPD@meta.data$celltype_sample))
colnames(cellnum) = c("group","cellnum")
rownames(cellnum) = cellnum$group
cellnum$group = NULL
pseudobulk = AddMetaData(pseudobulk, metadata = cellnum)

#filter cell numbers less than 10 in each celltype
pseudobulk = subset(pseudobulk, subset = cellnum >= 10)


