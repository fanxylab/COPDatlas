library(Seurat)
library(CellChat)
library(ComplexHeatmap)
library(tidyverse)

COPD = readRDS("./processed_data_R/COPD_Seurat_NMF.rds")
COPD.list = list("NMF1" = subset(COPD, subset = NMF_group == "P1"),
                 "NMF2" = subset(COPD, subset = NMF_group == "P2"))

CellChatDB = CellChatDB.human
CellChatDB.use = CellChatDB

for(i in c(1:2)) {
    COPD.list[[i]] = createCellChat(object = COPD.list[[i]], group.by = "celltype")
    COPD.list[[i]]@DB = CellChatDB.use
    COPD.list[[i]] = subsetData(COPD.list[[i]])
    COPD.list[[i]] = identifyOverExpressedGenes(COPD.list[[i]])
    COPD.list[[i]] = identifyOverExpressedInteractions(COPD.list[[i]])
    COPD.list[[i]] = computeCommunProb(COPD.list[[i]],population.size = FALSE)
    COPD.list[[i]] = filterCommunication(COPD.list[[i]], min.cells = 10)
    COPD.list[[i]] = computeCommunProbPathway(COPD.list[[i]])
    COPD.list[[i]] = aggregateNet(COPD.list[[i]])
    COPD.list[[i]] = netAnalysis_computeCentrality(COPD.list[[i]])
}

cellchat = mergeCellChat(COPD.list, add.names = names(COPD.list), cell.prefix = TRUE)

netVisual_heatmap(cellchat, width = 3, height = 3, measure = "weight", 
                  row.show = c("AT1","AT2","TRB Secretory",
                               "PreTB Secretory","Basal","Ciliated","Goblet"), 
                  col.show = c("B","Plasma","Treg","CD4+ T",
                               "CD8+ T","NK","cDC1","cDC2","Migratory DC","pDC",
                               "Classical monocyte","Non-classical monocyte",
                               "Alveolar macrophage","Monocyte-derived macrophage",
                               "Interstitial macrophage","Neutrophil"))

pathways.show = c("CXCL") 
weight.max = getMaxWeight(airway.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(airway.list)) {
  netVisual_aggregate(airway.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(names(airway.list)[i], pathways.show), color.use = airway_col,
                      sources.use = c("Goblet"), targets.use = c("B","Plasma","Treg","CD4+ T","CD8+ T","NK","cDC1","cDC2",
                                                        "Migratory DC","pDC","Classical monocyte","Non-classical monocyte",
                                                        "Alveolar macrophage","Monocyte-derived macrophage",
                                                        "Interstitial macrophage","Neutrophil") )
}

pathways.show = c("MHC-I") 
weight.max = getMaxWeight(airway.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(airway.list)) {
  netVisual_aggregate(airway.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(names(airway.list)[i], pathways.show), color.use = airway_col,
                      sources.use = c("Goblet"), targets.use = c("B","Plasma","Treg","CD4+ T","CD8+ T","NK","cDC1","cDC2",
                                                        "Migratory DC","pDC","Classical monocyte","Non-classical monocyte",
                                                        "Alveolar macrophage","Monocyte-derived macrophage",
                                                        "Interstitial macrophage","Neutrophil") )
}

netVisual_bubble(cellchat_airway, comparison = c(1, 2), signaling = c("IL6","CXCL","MHC-I","MIF"), 
                        color.text = group_col,
                        sources.use = c("Goblet"),
                        targets.use = c("B","Plasma","Treg","CD4+ T","CD8+ T","NK","cDC1","cDC2",
                                                        "Migratory DC","pDC","Classical monocyte","Non-classical monocyte",
                                                        "Alveolar macrophage","Monocyte-derived macrophage",
                                                        "Interstitial macrophage","Neutrophil"),
                        angle.x = 90)
