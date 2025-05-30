library(Seurat)
library(CellChat)
library(ComplexHeatmap)
library(tidyverse)

#load Seurat object
COPD = readRDS("processed_data_R/COPD_celltype_umap_0808_v3.rds")
Idents(COPD) = COPD@meta.data$Group
celltype = c("Immature AT1","Mature AT1","AT2","TRB Secretory","PreTB Secretory","Basal","Goblet","Ciliated","PNEC",
                 "Adventitial fibroblast","Alveolar fibroblast","Myofibroblast","Fibromyocyte","Pericyte","ASM",
                 "VSM","Mesothelial","Lymphatic","Aerocyte","gCap","Arterial","Venous","B","Plasma","Treg","CD4+ T",
                 "CD8+ T","Proliferating T","NKT","NK","Basophil","cDC1","cDC2","Migratory DC","pDC","Classical monocyte","Non-classical monocyte",
                 "Alveolar macrophage","Monocyte-derived macrophage","Interstitial macrophage","Mast","Neutrophil")
COPD@meta.data$celltype = factor(COPD@meta.data$celltype, levels = celltype)

COPD.list = list("HC" = subset(COPD, subset = Group == "HC"),
                 "HS" = subset(COPD, subset = Group == "HS"),
                 "Mild COPD" = subset(COPD, subset = Group == "Mild COPD"),
                 "Severe COPD" = subset(COPD, subset = Group == "Severe COPD"))

#select database
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

#Cellchat processing
for(i in c(1:4)) {
    COPD.list[[i]] <- createCellChat(object = COPD.list[[i]], group.by = "celltype")
    COPD.list[[i]]@DB <- CellChatDB.use
    COPD.list[[i]] <- subsetData(COPD.list[[i]])
    COPD.list[[i]] <- identifyOverExpressedGenes(COPD.list[[i]])
    COPD.list[[i]] <- identifyOverExpressedInteractions(COPD.list[[i]])
    COPD.list[[i]] <- computeCommunProb(COPD.list[[i]],population.size = FALSE)
    COPD.list[[i]] <- filterCommunication(COPD.list[[i]], min.cells = 10)
    COPD.list[[i]] <- computeCommunProbPathway(COPD.list[[i]])
    COPD.list[[i]] <- aggregateNet(COPD.list[[i]])
    COPD.list[[i]] <- netAnalysis_computeCentrality(COPD.list[[i]])
}

saveRDS(COPD.list, file = "./processed_data_R/cellchat_copd_list_group.rds")

#merge cellchat objects
cellchat = mergeCellChat(COPD.list, add.names = names(COPD.list), cell.prefix = TRUE)

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

#calculate cell-cell interaction numbers
num.link <- sapply(COPD.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(COPD.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(COPD.list[[i]], color.use = col, 
                                               title = names(COPD.list)[i], weight.MinMax = weight.MinMax) +
    xlim(0,80) + ylim(0,50)
}
patchwork::wrap_plots(plots = gg, ncol = 4)

#compare specific signaling pathways bewteen different groups
rankNet(cellchat, measure = "weight", mode = "comparison", stacked = T, show.raw = TRUE, 
               signaling = c("SELE","IGFBP","IL6","VCAM","FGF","WNT","EGF","TGFb","TNF","IFN-II","CXCL","CCL"),
               do.stat = FALSE, comparison =c(1,2,3,4), color.use = group_col, return.data = TRUE,do.flip = FALSE)

#draw L-R pairs in specific signaling pathways 
netVisual_bubble(cellchat, comparison = c(1,2,3,4), signaling = c("SELE","IGFBP","CXCL"), 
                 color.text = group_col, angle.x = 90,color.heatmap = "viridis",
                 sources.use = c("Venous","Arterial"),
                 targets.use = c("Classical monocyte","Non-classical monocyte",
                                "Alveolar macrophage","Monocyte-derived macrophage","CD4+ T","CD8+ T") )
