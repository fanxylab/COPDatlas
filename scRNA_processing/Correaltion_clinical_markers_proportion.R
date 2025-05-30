

celltype_ratio_df = read.csv("/datf/mazhuo/jupyter_notebook/COPD/tables/COPD_cellratio_new.csv", row.names = 1, check.names = FALSE)
newmeta = read.csv("tables/meta101.csv", row.names = 1, check.names = FALSE)
newmeta = merge(celltype_ratio_df, newmeta, by = "sample", all = TRUE)

cov_index = list()
index = c("Age","BMI","FEV1%pred","FEV1/FVC(%)","LAA","AWT-PI10")
for (x in index) {
    cov_index[[x]] = data.frame(group = x, celltype = celltype)
    for (cell in celltype)  {
    correlation = cor.test(newmeta[[x]], newmeta[[cell]], method = "pearson")
    p = correlation$p.value
    cor = correlation$estimate
    cov_index[[x]]$p_value[which(cov_index[[x]]$celltype == cell)] = p
    cov_index[[x]]$R[which(cov_index[[x]]$celltype == cell)] = cor
    }
}

cov$cellclass = "Immune"
cov$cellclass[which(cov$celltype %in% c("Immature AT1","Mature AT1","AT2","TRB Secretory","PreTB Secretory","Basal","Goblet","Ciliated","PNEC"))] = "Epithelial"
cov$cellclass[which(cov$celltype %in% c("Adventitial fibroblast","Alveolar fibroblast","Myofibroblast","Fibromyocyte","Pericyte","ASM","VSM"))] = "Mesenchymal"
cov$cellclass[which(cov$celltype %in% c("Lymphatic","Aerocyte","gCap","Arterial","Venous"))] = "Endothelial"
cov$cellclass = factor(cov$cellclass, levels = c("Epithelial","Mesenchymal","Endothelial","Immune"))
cov$celltype = factor(cov$celltype, levels = c("Immature AT1","Mature AT1","AT2","TRB Secretory","PreTB Secretory","Basal","Goblet","Ciliated","PNEC",
                                               "Adventitial fibroblast","Alveolar fibroblast","Myofibroblast","Fibromyocyte","Pericyte","ASM","VSM",
                                               "Lymphatic","Aerocyte","gCap","Arterial","Venous","B","Plasma","Treg","CD4+ T","CD8+ T","Proliferating T","NKT","NK","cDC1","cDC2","Migratory DC",
                                               "pDC","Classical monocyte","Non-classical monocyte","Alveolar macrophage","Monocyte-derived macrophage","Interstitial macrophage",
                                               "Mast","Neutrophil","Basophil"))
cov$group = factor(cov$group, levels = c("AWT-PI10","LAA","FEV1/FVC(%)","FEV1%pred","BMI","Age"))

ggplot(data = cov, aes(x = celltype, y = group, size = `-log10(p_value)`, col = R)) +
 geom_point() + theme_bw() + xlab("") + ylab("") +  
 facet_grid(. ~ cellclass, scales = "free_x", space = "free_x") +
 scale_color_gradient2(high="#E31A1C",mid = "lightgrey",low ="#5A8FCA", midpoint = 0) +
 theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
       axis.text.y = element_text(color = "black"),       
       panel.grid.major = element_blank(), panel.grid.minor =element_blank())
