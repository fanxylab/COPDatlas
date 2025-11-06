library(tidyverse)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)
library(gprofiler2)
library(GO.db)
library(AnnotationDbi)
library(reshape2)
library(ComplexHeatmap)
library(colorRamp2)

#Fig. 5d
gene_cor_FEV1_all = read.csv(file = "./tables/gene_cor_FEV1_all_241219.csv", row.names = 1)
celltype = unique(gene_cor_FEV1_all$celltype)

term = c("GO:0006119","GO:0043484","GO:0043484","GO:0016055","GO:0097193","GO:0019915",
         "GO:0097193","GO:0002757","GO:0038202","GO:0007015","GO:0019882","GO:0010506",
         "GO:0043410","GO:0070555","GO:0006979","GO:0007219","GO:0030509","GO:0007179")
term_list = list()
for (x in term) {
 term_list[x] = Term(x)
    }
term_all = do.call(rbind, term_list)
term_all = as.data.frame(term_all)
colnames(term_all) = "Term"
term_all$GOALL = rownames(term_all)

GO_genes = bitr(geneID = term,
                        fromType = "GOALL",
                        toType = "SYMBOL",
                        OrgDb = org.Hs.eg.db)
GO_genes = left_join(GO_genes, term_all, by = "GOALL")

gse_all = list()
for(x in celltype) {
    gene = gene_cor_FEV1_all[which(gene_cor_FEV1_all$celltype == x),]
    gene_list = gene$R
    names(gene_list) = gene$gene
    gene_list = sort(gene_list, decreasing = TRUE)
    gseaResult = GSEA(
       geneList = gene_list,
       exponent = 1,
       minGSSize = 10,
       maxGSSize = 500,
       pvalueCutoff = 1, 
       TERM2GENE = data.frame(term = GO_genes$Term,
                           gene = GO_genes$SYMBOL),
       seed = TRUE
      )
    gse = as.data.frame(gseaResult@result)
    gse$celltype = x
    gse_all[[x]] = gse
    }
gse_df = do.call(rbind, gse_all)

celltype = c('Immature AT1','AT1','AT2','TRB Secretory','PreTB Secretory','Basal','Goblet','Ciliated',
                    'Adventitial fibroblast','Alveolar fibroblast','Pericyte','ASM',
                    'Lymphatic','Aerocyte','gCap','Arterial','Venous','B','Plasma','Treg','CD4+ T','CD8+ T','NKT',
                    'NK','cDC2','pDC','Classical monocyte','Non-classical monocyte','Alveolar macrophage',
                    'Monocyte-derived macrophage','Interstitial macrophage','Mast','Neutrophil')

gse_df$celltype = factor(gse_df$celltype, levels = celltype)
gse_matrix = dcast(gse_df, Description ~ celltype, value.var = "NES")
rownames(gse_matrix) = gse_matrix$Description
gse_matrix$Description = NULL
gse_matrix[is.na(gse_matrix)] = 0

fdr_matrix = dcast(gse_df, Description ~ celltype, value.var = "pvalue")
rownames(fdr_matrix) = fdr_matrix$Description
fdr_matrix$Description = NULL
fdr_matrix[is.na(fdr_matrix)] = 1

significant_label = ifelse(fdr_matrix < 0.01, "**",
                           ifelse(fdr_matrix < 0.05, "*", ""))

Heatmap(gse_matrix, name = "NES", 
        cluster_rows = FALSE, cluster_columns = FALSE,
        border = TRUE, row_names_side = c("left"),
        col = colorRamp2(c(-3, 0, 3), c('#5A8FCA','#F2F2F0','#E31A1C')),
        cell_fun = function(j, i, x, y, width, height, fill) {
             grid.text(significant_label[i,j], x, y, gp = gpar(fontsize = 10))
} )
