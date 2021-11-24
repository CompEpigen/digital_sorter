library(Seurat)
#library(SeuratDisk)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)

wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)
ob = readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/song_2019.rds")
normal = read.csv("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/procdata/maker_gene_expression_in_normal_lung.csv")
#ob@meta.data$disease[!(ob@meta.data$disease=="adjacent normal"|ob@meta.data$disease=="normal"|ob@meta.data$disease=="nsclc")]="nsclc_meta"
table(ob@meta.data$disease)
#adjacent normal 6598 #nsclc 4883
table(ob@meta.data$donor)
#P1_Normal  P1_Tumor P2_Normal  P2_Tumor P3_Normal  P3_Tumor P4_Normal  P4_Tumor 
#  2495      1832      1409      1300       667       328      2027      1423 


cols = c("steelblue","darkred","gold","coral2")

# CD45 (PTPRC) ####
features=c("PTPRC")
m1 =VlnPlot(ob, features, split.by = "disease", group.by = "annotation.l2", cols=cols,
            sort = TRUE,pt.size = 0, combine = FALSE)
m1

d = m1[[1]]$data
mu = aggregate(x=d$PTPRC,
               by=list(d$ident,d$split),
               FUN=mean)
summary(mu)
normal_PTPRC <- mean(normal$PTPRC, na.rm = T)

set1 = unique(sort(normal[normal$PTPRC >=normal_PTPRC,]$Cell.type))
set2 = levels(set1)[!levels(set1) %in% set1]
meta = ob$annotation.l2
meta[meta %in% set1] = paste0(features, "+")
meta[meta %in% set2] = paste0(features, "-")


ob <- AddMetaData(object = ob, metadata = meta, col.name = "split1.1")


# EPCAM ####
ob_split <- SplitObject(ob, split.by = "split1.1")
ob_split_PTPRC0 = ob_split$`PTPRC-`

features=c("EPCAM")
m1 =VlnPlot(ob_split_PTPRC0, features, split.by = "disease", group.by = "annotation.l2", 
            cols = cols,sort = TRUE,pt.size = 0, combine = FALSE)
m1
d = m1[[1]]$data
mu = aggregate(x=d$EPCAM,
               by=list(d$ident,d$split),
               FUN=mean, na.action=)
summary(mu)
normal_EPCAM <-  mean(normal$EPCAM)

set1 = unique(sort(normal[normal$EPCAM >=normal_EPCAM,]$Cell.type))
set2 = levels(set1)[!levels(set1) %in% set1]
meta = ob$annotation.l2
meta[meta %in% set1] = paste0(features, "+")
meta[meta %in% set2] = paste0(features, "-")

ob_split_PTPRC0 <- AddMetaData(object = ob_split_PTPRC0, metadata = meta, col.name = "split2")
ob2 <- AddMetaData(object = ob, metadata = meta, col.name = "split2")
ob_split = SplitObject(ob, split.by = "disease")



## Level 3
ob3 <- SplitObject(ob2, split.by = "split2")
ob_split_EPCAM0 = ob3$`EPCAM-`

features=c("PECAM1")
m1 =VlnPlot(ob_split_EPCAM0, features, split.by = "disease", group.by = "annotation.l2", 
            cols = cols,sort = TRUE,pt.size = 0, combine = FALSE)

m1
d = m1[[1]]$data
mu = aggregate(x=d$PECAM1,
               by=list(d$ident,d$split),
               FUN=mean, na.action=)
summary(mu)
normal_PECAM1 <- mean(normal$PECAM1)

set1 = unique(sort(normal[normal$PECAM1 >=normal_PECAM1,]$Cell.type))
set2 = levels(set1)[!levels(set1) %in% set1]
meta = ob$annotation.l2
meta[meta %in% set1] = paste0(features, "+")
meta[meta %in% set2] = paste0(features, "-")

ob3 <- AddMetaData(object = ob2, metadata = meta, col.name = "split3")

saveRDS(ob3,"/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/song_2019_addsplits.rds")


####### CL ##############
CL <- read.csv("procdata/CL_cell_names_lee.csv")
CL <- CL[,3]
CL <- unique(sort(CL))
anno1 <- ob@meta.data[["annotation.l1"]]
anno1 <- unique(sort(anno1))

anno2 <- ob@meta.data[["annotation.l2"]]
anno2 <- unique(sort(anno2))

n <- max(length(anno1), length(anno2),length(CL))
length(anno1) <- n                      
length(anno2) <- n

cell_anno <- cbind(CL,anno1,anno2)
colnames(cell_anno) <- c("Cell ontology","annotation.l1","annotation.l2")
write.csv(cell_anno,"dirty_scripts/compare_cell_anno.csv")
