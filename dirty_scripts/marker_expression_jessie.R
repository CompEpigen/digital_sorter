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
#ob = readRDS("/omics/groups/OE0219/internal/MJMC/P01_NSCLC/P01.2_enriched_cell_components/analysis/azimuth-meta-analysis/human_lung/cellxgene/lung_mapped_cellxgene.rds")

## deal with normal table
#normal = read.csv("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/procdata/maker_gene_expression_in_normal_lung.csv")
#normal_EPCAM <- readRDS("dirty_scripts/normal_mean_EPCAM.rds")
#normal_EPCAM <- normal_EPCAM[,c(1,3)]
#normal2 <- merge(normal, normal_EPCAM, by.x = "Cell.type", by.y = "Group.1", all.x =T)
#normal_PECAM1 <- readRDS("dirty_scripts/normal_mean_PECAM1.rds")
#normal_PECAM1 <- normal_PECAM1[,c(1,3)]
#normal3 <- merge(normal2, normal_PECAM1, by.x = "Cell.type", by.y = "Group.1", all.x =T)
#normal3 <- normal3[,c(1,3,4,7,8)]
#colnames(normal3) <- c("Cell.type","source","PTPRC","EPCAM","PECAM1")
#saveRDS(normal3,"procdata/maker_gene_expression_in_normal_lung_jessie.RDS")
normal <- readRDS("procdata/maker_gene_expression_in_normal_lung_jessie.RDS")

cols = c("steelblue","darkred","gold","coral2")


## CD45 (PTPRC) ####
features=c("PTPRC")
m1 =VlnPlot(ob, features, split.by = "disease", group.by = "annotation.l2", cols=cols,
            sort = TRUE,pt.size = 0, combine = FALSE)
m1

d = m1[[1]]$data
cell <- as.character(unique(sort(ob@meta.data[["annotation.l2"]])))
normal_sub <- normal[normal$Cell.type %in% cell,]
normal_PTPRC <- mean(normal_sub$PTPRC, na.rm = T)

PTPRC_group <- c()
for(i in 1:length(normal$PTPRC)){
  if(normal$PTPRC[i] >=normal_PTPRC){PTPRC_group <- append(PTPRC_group,"PTPRC+")}
  else{PTPRC_group <- append(PTPRC_group,"PTPRC-")}
}

normal$PTPRC_group <- PTPRC_group
set1 <- unique(sort(normal[normal$PTPRC_group =="PTPRC+",]$Cell.type))
set2 <- unique(sort(normal[normal$PTPRC_group =="PTPRC-",]$Cell.type))
meta = ob$annotation.l2
meta[meta %in% set1] = paste0(features, "+")
meta[meta %in% set2] = paste0(features, "-")

ob <- AddMetaData(object = ob, metadata = meta, col.name = "split1.1")


## EPCAM ####
ob2 <- SplitObject(ob, split.by = "split1.1")
ob2 = ob2$`PTPRC-`

features=c("EPCAM")
m1 =VlnPlot(ob2, features, split.by = "disease", group.by = "annotation.l2", 
            cols = cols,sort = TRUE,pt.size = 0, combine = FALSE)
m1
d = m1[[1]]$data
cell <- as.character(unique(sort(ob2@meta.data[["annotation.l2"]])))
normal_sub <- normal[normal$Cell.type %in% cell,]
normal_EPCAM <-  mean(normal_sub$EPCAM) #median not mean #use mean: club cell would be EPCAM-
EPCAM_group <- c()
for(i in 1:length(normal$EPCAM)){
  if (is.na(normal$EPCAM[i])) {EPCAM_group <- append(EPCAM_group,"NA")}
  else if(normal$EPCAM[i] >=normal_EPCAM){EPCAM_group <- append(EPCAM_group,"EPCAM+")}
  else{EPCAM_group <- append(EPCAM_group,"EPCAM-")}
}


normal$EPCAM_group <- EPCAM_group
set1 <- unique(sort(normal[normal$EPCAM_group =="EPCAM+",]$Cell.type))
set2 <- unique(sort(normal[normal$EPCAM_group =="EPCAM-",]$Cell.type))
meta = ob$annotation.l2
meta[meta %in% set1] = paste0(features, "+")
meta[meta %in% set2] = paste0(features, "-")

ob <- AddMetaData(object = ob, metadata = meta, col.name = "split2")
ob2 <- AddMetaData(object = ob2, metadata = meta, col.name = "split2")




## PECAM1 ####
ob3 <- SplitObject(ob2, split.by = "split2")
ob3 = ob3$`EPCAM-`

features=c("PECAM1")
m1 =VlnPlot(ob3, features, split.by = "disease", group.by = "annotation.l2", 
            cols = cols,sort = TRUE,pt.size = 0, combine = FALSE)

m1
d = m1[[1]]$data
cell <- as.character(unique(sort(d$ident)))
normal_sub <- normal[normal$Cell.type %in% cell,]
normal_PECAM1 <-  median(normal_sub$PECAM1) #median instead of mean
PECAM1_group <- c()
for(i in 1:length(normal$PECAM1)){
  if (is.na(normal$PECAM1[i])) {PECAM1_group <- append(PECAM1_group,"NA")}
  else if(normal$PECAM1[i] >=normal_PECAM1){PECAM1_group <- append(PECAM1_group,"PECAM1+")}
  else{PECAM1_group <- append(PECAM1_group,"PECAM1-")}
}


normal$PECAM1_group <- PECAM1_group
set1 <- unique(sort(normal[normal$PECAM1_group =="PECAM1+",]$Cell.type))
set2 <- unique(sort(normal[normal$PECAM1_group =="PECAM1-",]$Cell.type))
meta = ob$annotation.l2
meta[meta %in% set1] = paste0(features, "+")
meta[meta %in% set2] = paste0(features, "-")

ob <- AddMetaData(object = ob, metadata = meta, col.name = "split3")

saveRDS(ob,"/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/song_2019_addsplits.rds")

normal$X <- NULL

write.csv(normal,"procdata/maker_gene_expression_in_normal_lung_add_groups.csv")

