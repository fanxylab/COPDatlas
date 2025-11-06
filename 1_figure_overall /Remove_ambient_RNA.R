library(SoupX)
library(DropletUtils)

##Load data from cellranger output and remove ambient RNA contamination
input_dir <- "/datf/mazhuo/cellranger/"

for (id in c("COPD2","COPD3","COPD4","COPD5","COPD7","COPD8","COPD9","COPD10","COPD11","COPD12","COPD15","COPD16",
             "COPD21","COPD23","COPD26","COPD27","COPD28","COPD29","COPD30","COPD31","COPD32","COPD38","COPD39",
             "COPD40","COPD41","COPD45","COPD46","COPD50","COPD52","COPD53","COPD54","COPD56","COPD57","COPD58",
             "COPD64","COPD71","COPD72","COPD73","COPD77","COPD78","COPD79","COPD81","COPD82","COPD83",
             "COPD87","COPD88","COPD90","COPD91","COPD92","COPD93","COPD94","COPD95","COPD97","COPD99",
             "COPD101","COPD102","COPD103","COPD104","COPD107","COPD128",
             "HC2","HC3","HC5","HC6","HC10","HC24","HC26","HC27","HC29","HC30","HC31","HC32","HC36","HC38","HC40",
             "HC46","HC47","HC49","HC50","HS1","HS2","HS6","HS7","HS14","HS21","HS22","HS24","HS25","HS26","HS27",
             "HS28","HS29","HS30","HS31","HS32","HS34","HS35","HS41","HS43","HS44","HS46")) 
{
  data <- load10X(paste0(input_dir,id,"/outs/"))
  data <- autoEstCont(data)
  data <- adjustCounts(data, roundToInt = T)
  write10xCounts(paste0("/datf/mazhuo/data/COPD/Matrix_after_SoupX/",id), data, version="3")
}
