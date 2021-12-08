## inspect the merge seurat object ####
if(T){
  library(Seurat)
  #library(SeuratDisk)
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  
  wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
  setwd(wd)
  dataset <- readRDS("rawdata/lung_mapped_cellxgene_fixed_addsplits.rds")
  split <- as.character(unique(sort(dataset@meta.data$dataset_origin)))
  ob_all <- SplitObject(dataset, split.by = "dataset_origin")
  saveRDS(ob_all,"rawdata/lung_mapped_cellxgene_seuratsplit_fixed_7datasets.rds")
  
  ob1 <- ob_all[[1]]
  ob2 <- ob_all[[2]]
  ob3 <- ob_all[[3]]
  ob4 <- ob_all[[4]]
  ob5 <- ob_all[[5]]
  ob6 <- ob_all[[6]]
  ob7 <- ob_all[[7]]
  split <- names(ob_all)
  list_samples_disease <- list(
    as.character(unique(sort(ob1@meta.data[["disease"]]))),
    as.character(unique(sort(ob2@meta.data[["disease"]]))),
    as.character(unique(sort(ob3@meta.data[["disease"]]))),
    as.character(unique(sort(ob4@meta.data[["disease"]]))),
    as.character(unique(sort(ob5@meta.data[["disease"]]))),
    as.character(unique(sort(ob6@meta.data[["disease"]]))),
    as.character(unique(sort(ob7@meta.data[["disease"]])))
  )
  names(list_samples_disease) <-names(ob_all)
  
  saveRDS(list_samples_disease,"rawdata/LIST_lung_seuratsplit_7datasets_samples_disease.rds")
}


if(F){ 
## Add splits into the merged seurat object ####

wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)

ob = readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/lung_mapped_cellxgene.rds")

## deal with normal table
if(F){ 
  normal = read.csv("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/procdata/maker_gene_expression_in_normal_lung_fixed.csv")
  normal_EPCAM <- readRDS("dirty_scripts/normal_mean_EPCAM.rds")
  normal_EPCAM <- normal_EPCAM[,c(1,3)]
  normal2 <- merge(normal, normal_EPCAM, by.x = "Cell.type", by.y = "Group.1", all.x =T)
  normal_PECAM1 <- readRDS("dirty_scripts/normal_mean_PECAM1.rds")
  normal_PECAM1 <- normal_PECAM1[,c(1,3)]
  normal3 <- merge(normal2, normal_PECAM1, by.x = "Cell.type", by.y = "Group.1", all.x =T)
  normal3 <- normal3[,c(1,3,4,7,8)]
  colnames(normal3) <- c("Cell.type","source","PTPRC","EPCAM","PECAM1")
  normal_MME <- readRDS("dirty_scripts/normal_mean_MME.rds")
  normal_MME <- normal_MME[,c(1,3)]
  normal4 <- merge(normal3, normal_MME, by.x = "Cell.type", by.y = "Group.1", all.x =T)
  colnames(normal4) <- c("Cell.type","source","PTPRC","EPCAM","PECAM1","MME")
  saveRDS(normal4,"procdata/maker_gene_expression_in_normal_lung_jessie_fixed.RDS")
}
#defined with the normal tissues from dataset travaglini_2020
normal <- readRDS("procdata/maker_gene_expression_in_normal_lung_jessie_fixed.RDS")
diseases <- as.character(unique(sort(ob@meta.data$disease)))
cells <- as.character(unique(sort(ob@meta.data$annotation.l2)))
cols = c("steelblue","darkred","gold","coral2","pink","yellowgreen","violet","blue","grey","black")


## CD45 (PTPRC) ####
features=c("PTPRC")
m1 =VlnPlot(ob, features, split.by = "disease", group.by = "annotation.l2", cols=cols,
            sort = TRUE,pt.size = 0, combine = FALSE)
#m1

d = m1[[1]]$data

normal_PTPRC <- mean(normal$PTPRC, na.rm = T)

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

ob <- AddMetaData(object = ob, metadata = meta, col.name = "split1")


## EPCAM ####
ob2 <- SplitObject(ob, split.by = "split1")
ob2 = ob2$`PTPRC-`

features=c("EPCAM")
m1 =VlnPlot(ob2, features, split.by = "disease", group.by = "annotation.l2", 
            cols = cols,sort = TRUE,pt.size = 0, combine = FALSE)
#m1
d = m1[[1]]$data

normal_sub <- normal[normal$PTPRC_group== "PTPRC-",]
normal_EPCAM <-  mean(normal_sub$EPCAM,na.rm=T) 
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

#m1
d = m1[[1]]$data

normal_sub <- normal[normal$EPCAM_group== "EPCAM-",]
normal_PECAM1 <-  mean(normal_sub$PECAM1,na.rm=T)  
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
ob3 <- AddMetaData(object = ob3, metadata = meta, col.name = "split3")


## MME ####
ob4 <- SplitObject(ob3, split.by = "split3")
ob4 = ob4$`PECAM1-`

features=c("MME")
m1 =VlnPlot(ob4, features, split.by = "disease", group.by = "annotation.l2", 
            cols = cols,sort = TRUE,pt.size = 0, combine = FALSE)

#m1
d = m1[[1]]$data
normal_sub <- normal[normal$PECAM1_group== "PECAM1-",]
normal_MME <-  mean(normal_sub$MME,na.rm=T) 
MME_group <- c()
for(i in 1:length(normal$MME)){
  if (is.na(normal$MME[i])) {MME_group <- append(MME_group,"NA")}
  else if(normal$MME[i] >=normal_MME){MME_group <- append(MME_group,"MME+")}
  else{MME_group <- append(MME_group,"MME-")}
}


normal$MME_group <- MME_group
set1 <- unique(sort(normal[normal$MME_group =="MME+",]$Cell.type))
set2 <- unique(sort(normal[normal$MME_group =="MME-",]$Cell.type))
meta = ob$annotation.l2
meta[meta %in% set1] = paste0(features, "+")
meta[meta %in% set2] = paste0(features, "-")

ob <- AddMetaData(object = ob, metadata = meta, col.name = "split4")



saveRDS(ob,"/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/lung_mapped_cellxgene_fixed_addsplits.rds")
}
