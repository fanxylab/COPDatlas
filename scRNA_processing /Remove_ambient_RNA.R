library(SoupX)
library(DropletUtils)

##Load data from cellranger output and remove ambient RNA contamination
input_dir <- "/datf/mazhuo/cellranger/"

for (id in c("COPD1","COPD2","COPD3","COPD4","COPD5","COPD6","COPD7","COPD8","COPD9","COPD10","COPD11","COPD12",
              "COPD15","COPD16","COPD21","COPD23","COPD24","COPD26","COPD27","COPD28","COPD29","COPD30",
              "COPD31","COPD32","COPD33","COPD38","COPD39","COPD40","COPD41","COPD45","COPD46","COPD50",
              "COPD52","COPD53","COPD54","COPD56","COPD57","COPD58",
              "HS1","HS2","HS6","HS7","PC1","PC2","HC1","HC2","HC3","HC5","HC6","HC9","HC10","JK03","JK04")) 
{
  data <- load10X(paste0(input_dir,id,"/outs/"))
  data <- autoEstCont(data)
  data <- adjustCounts(data, roundToInt = T)
  write10xCounts(paste0("/datf/mazhuo/data/COPD/Matrix_after_SoupX/",id), data, version="3")
}
