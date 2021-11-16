wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/dirty_scripts"
setwd(wd)

library(Seurat)
#library(SeuratDisk)
library(dplyr)
library(cerebroApp)
library(grid)
library(gridExtra)
library(ggplot2)

ob = readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/song_2019_addsplits.rds")
cols = c("steelblue","darkred","gold","coral2")
ob_split <- SplitObject(ob, split.by = "disease")
ob1 <- SplitObject(ob_split[[1]], split.by = "split1")
features=c("CD8A","CD8B","CD4", "MRC1","CD14","SIGLEC1")
DotPlot(ob1$`PTPRC+`, features = features, group.by = "annotation.l2") + RotatedAxis()



ob2 <- SplitObject(ob_split[[1]], split.by = "split1.1")
features=c("CD8A","CD8B","CD4", "MRC1","CD14","SIGLEC1")
DotPlot(ob2$`PTPRC+`, features = features, group.by = "annotation.l2") + RotatedAxis()
