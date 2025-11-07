library(Seurat)
library(monocle3)
library(tidydr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

COPD = readRDS("./processed_data_R/COPD_celltype.rds")
preTB_meta = read.csv("./tables/preTB_meta.csv", row.names = 1, check.names = FALSE)
preTB = subset(COPD, subset = celltype %in% c("AT1","Immature AT1","TRB Secretory","AT2",
                                            "PreTB Secretory") )
preTB@meta.data = preTB_meta
UMAP = as.matrix(preTB@meta.data[,c("UMAP1","UMAP2")])
colnames(UMAP) = c("umap_1","umap_2")
preTB@reductions$umap@cell.embeddings = UMAP

data = preTB@assays$RNA@counts
cell_metadata = preTB@meta.data
gene_annotation = data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) = rownames(data)
cds = new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds = preprocess_cds(cds, num_dim = 50)
cds = reduce_dimension(cds, preprocess_method = "PCA")
cds.embed = cds@int_colData$reducedDims$UMAP
int.embed = Embeddings(preTB, reduction = "umap")
int.embed = int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP = int.embed

cds = cluster_cells(cds)
cds = learn_graph(cds,use_partition = F, learn_graph_control = list(minimal_branch_len = 2))

UMAP = as.data.frame(UMAP)
UMAP = UMAP[order(-UMAP$umap_2),]
cds = order_cells(cds, root_cells = rownames(UMAP[1,]))

plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_roots = F,cell_size =0.5,
                label_leaves = F,  label_branch_points = F,trajectory_graph_color = "grey", rasterize = TRUE,
                trajectory_graph_segment_size = 0.4) +
  theme_dr() + theme(panel.grid = element_blank()) +
  scale_color_viridis_c() + labs(color = "Pseudotime")

celltype_col = setNames(c("#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#b3bd22", "#17becf", "#aec7e8"),
                   c("Immature AT1","AT1","AT2","TRB Secretory","PreTB Secretory","Basal","Goblet","Ciliated","PNEC"))

plot_cells(cds, color_cells_by = "celltype", label_cell_groups = F, label_roots = F,cell_size =0.5,
           label_leaves = F,  label_branch_points = F,trajectory_graph_color = "grey", rasterize = TRUE,
                trajectory_graph_segment_size = 0.4) +
      scale_color_manual(values = celltype_col) +
      theme_dr() + theme(panel.grid = element_blank()) 

pseudotime = pseudotime(cds) %>% as.data.frame()
colnames(pseudotime)[1] = "pseudotime"
preTB = AddMetaData(preTB, metadata = pseudotime)

group_col = setNames(c("#51a246","#3194cc","#794191","#684b3b"),
                     c("HC","HS","Mild COPD","Severe COPD") )

preTB_sub = preTB@meta.data[which(preTB@meta.data$pseudotime < 15 &
                                  preTB@meta.data$celltype %in% c("PreTB Secretory","TRB Secretory") ),]

ggplot(preTB_sub,aes(pseudotime,fill=celltype, color=celltype)) +  
    geom_density(alpha = 0.5,size=1.2)  +
    scale_fill_manual(values=celltype_col)  + 
    scale_color_manual(values=celltype_col) +
    labs(x = "Pseudotime", y = "Density") + 
    theme_classic() +
    scale_y_continuous(expand = c(0,0)) + 
    theme(axis.text=element_text(colour='black',size=10))
