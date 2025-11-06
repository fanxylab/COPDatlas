library(networkD3)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggalluvial)

celltypist = read.csv("./tables/COPD_meta_celltypist.csv", row.names = 1, check.names = FALSE)

celltype = c("Immature AT1","AT1","AT2","TRB Secretory","PreTB Secretory","Basal","Goblet","Ciliated","PNEC",
"Adventitial fibroblast","Alveolar fibroblast","Myofibroblast","Fibromyocyte","Pericyte","ASM",
"VSM","Mesothelial","Lymphatic","Aerocyte","gCap","Arterial","Venous","B","Plasma","Treg","CD4+ T",
"CD8+ T","Proliferating T","NKT","NK","cDC1","cDC2","Migratory DC","pDC","Classical monocyte","Non-classical monocyte",
"Alveolar macrophage","Monocyte-derived macrophage","Interstitial macrophage","Mast","Neutrophil","Basophil")

celltype_col = setNames(c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896","#9467bd", "#c5b0d5", 
                     "#8c564b", "#c49c94", "#e377c2", "#f786d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d","#17becf","#9edae5", 
                     "#6baed6", "#fd8d3c", "#74c476", "#9e9ac7", "#969696", "#5254a3", "#8ca252", "#bd9e39", "#ad494a", "#a55194",
                     "#FF410D", "#6EE2FF", "#F7C530", "#95cc5e", "#d0dfe6", "#f79d1e", "#748aa6", "#cc0c00", "#5c88da", "#84bd00",
                     "#ffcd00", "#7c878e"),
                c('PNEC','Basal','TRB Secretory','AT2','PreTB Secretory','Goblet','Ciliated',
                  'Immature AT1','AT1','Mesothelial','VSM','ASM','Pericyte','Myofibroblast','Fibromyocyte',
                  'Alveolar fibroblast','Adventitial fibroblast','Lymphatic','gCap','Aerocyte','Arterial',
                  'Venous','Migratory DC','cDC1','pDC','Basophil','Mast','Plasma','B',
                  'Interstitial macrophage','Monocyte-derived macrophage','Alveolar macrophage','cDC2',
                  'Classical monocyte','Non-classical monocyte','Neutrophil','CD4+ T','Treg','CD8+ T',
                  'Proliferating T','NKT','NK'))

type = c("AT1","AT2","AT0","pre-TB secretory","Basal resting","Multiciliated (non-nasal)",
         "Neuroendocrine","Adventitial fibroblasts","Alveolar fibroblasts","Peribronchial fibroblasts",
         "Myofibroblasts","Pericytes","Smooth muscle","Mesothelium","Lymphatic EC differentiating",
         "EC aerocyte capillary","EC arterial","EC general capillary","EC venous systemic",
         "EC venous pulmonary","B cells","Plasma cells","CD4 T cells","CD8 T cells","NK cells",
         "DC1","DC2","Migratory DCs","Plasmacytoid DCs","Classical monocytes","Non-classical monocytes",
         "Alveolar macrophages","Alveolar Mph proliferating","Monocyte-derived Mph","Interstitial Mph perivascular","Mast cells")

celltypist$celltype = factor(celltypist$celltype, levels = celltype)
celltypist$majority_voting = factor(celltypist$majority_voting, levels = type)

sankey_data = sankey_data %>%
  group_by(celltype, majority_voting) %>%
  mutate(Count = 1) %>%  
  ungroup()

ggplot(sankey_data,
       aes(axis1 = celltype,        
           axis2 = majority_voting, 
           y = Count)) +            
  geom_alluvium(aes(fill = celltype),  
                alpha = 0.7,           
                width = 0.3) +        
  geom_stratum(width = 0.3,          
               fill = "grey90",       
               color = "black") +    
  geom_text(stat = "stratum",        
            aes(label = after_stat(stratum)),
            size = 3) +
  scale_x_discrete(limits = c("celltype", "majority_voting")) +  
  scale_y_continuous(breaks = NULL) +  
  theme_minimal() + scale_fill_manual(values = celltype_col) +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) 
