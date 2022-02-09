library(Seurat)
#library(SeuratDisk)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)  
print("set wd")

if(F){ 
## download one of the dataset for testing SurfaceGenie  ####
ob1 <- ob_all[[1]]
ob_split <- SplitObject(ob1, split.by = "disease")
ob_split2 = ob_split[[1]]
count <- as.data.frame(ob_split2@assays$RNA@counts)
count[1:3,1:3]
count$Accession <- row.names(count)
count <- count %>% relocate(Accession)
count[1:3,1:3]
count$Accession <- gsub('\\..*','',count$Accession)
count[1:3,1:3]
row.names(count) <- NULL
write.csv(count,"SurfaceGenie/song_2019_adjNormal_countTable_for_SurfaceGenie.csv")
 }
## download one of the dataset for testing SurfaceGenie  ####
SPC <- read.delim("SurfaceGenie/ref_SPC_by_Source_sprot.txt", sep=",") #from https://github.com/GundryLab/SurfaceGenie/tree/master/ref
SPC <- SPC %>% filter(SPC >= 3)
SPC2 <- SPC[,1:2]


library("UniprotR")
print("UniprotR")
ids <- SPC2$Accession
library(glue)
ids2 <-glue_collapse(ids, sep = " ")
#paste ids2 to https://www.uniprot.org/uploadlists/ §20220110
#download mapping table from the website and save as SurfaceGenie/SPC_highScore_geneID.txt
SPC_high <- read.delim("SurfaceGenie/SPC_highScore_geneID.txt") 
cell.surface.marker <- readRDS("digital_sorter/data/cell.surface.marker.rds") #883
new_cell_marker_gene <-SPC_high[,2] #2229
union_marker <- union(new_cell_marker_gene,cell.surface.marker) #2557
union_marker = union_marker[!grepl("HLA",union_marker)] #2536
saveRDS(union_marker,"SurfaceGenie/Union.cell.surface.marker.rds")
saveRDS(new_cell_marker_gene,"SurfaceGenie/SurfaceGenie.cell.surface.marker.rds")
