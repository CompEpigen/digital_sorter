# Workflow of adding splits AND marker selection of our own lung patient data ####
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)

# Set working directory ####
wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)
rawdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/patient_hlcma"

outdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/patient_hlcma_addsplits"
source("dirty_scripts/getMarkers.R")
# Load functions
GetfileNames <- function(fileDir, pattern = ".rds"){
  filenames <- list.files(fileDir, pattern= pattern, full.names=TRUE)
  filename <- sub(paste0(".*",fileDir,"/"), "", filenames)
  filename <- sub(pattern, "", filename)  
  return(filename)
}


# Read files ####
#ReadRDSFiles(rawdir, envir = .GlobalEnv)
filename <- GetfileNames(rawdir, pattern = ".rds")

## cell types grouped by master markers' expression created in 01_marker_expression_jessie.R
normal <- read.csv("procdata/maker_gene_expression_in_normal_lung_add_groups_fixed.csv")

cell.surface.marker <-readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/SurfaceGenie/Union.cell.surface.marker.rds")



if(F){
  ob <- readRDS(paste0(rawdir,"/",filename[1],".rds"))
  sampledata <- ob@meta.data
  dataset <- unique(sort(sampledata$dataset_origin))
  sampletype <- unique(sort(sampledata$disease))
  
  ob <- SplitObject(ob, split.by = "dataset_origin")
  ob <- ob[["hlcma0001"]]
  # add split
  cols = c("steelblue","darkred","gold","coral2","pink","yellowgreen","violet","blue","grey","black")
  
  
  ## CD45 (PTPRC) ####
  features=c("PTPRC")
  m1 =VlnPlot(ob, features, split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
  m1
  
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
  
  
  
  saveRDS(ob,"/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/patient_hlcma/patient_hlcma_addsplits.rds")
  
}



ob <-readRDS(paste0(rawdir,"/patient_hlcma_addsplits.rds"))

###Compare markers between/among datasets####
split <- as.character(unique(sort(ob@meta.data$disease)))
ob_split <- SplitObject(ob, split.by = "disease")


#k <- 1L
for(k in 1:length(split)){ 

print(split[k])
ob_split2 = ob_split[[split[k]]]
ob_split3 <- SplitObject(ob_split2, split.by = "split1")
ob_split3 = ob_split3$`PTPRC-`
ob_split4 <- SplitObject(ob_split3, split.by = "split2")
ob_split4 = ob_split4$`EPCAM+`

seurat <- getMarkers(
  ob_split2,assay = 'RNA',
  organism = 'hg',
  groups = c('annotation.l2'),
  name = 'cerebro_seurat',
  only_pos = F,
  min_pct = 0.5,
  thresh_logFC = 1,
  thresh_p_val = 0.05,
  test = 'wilcox', 
  verbose = TRUE
)
marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]] #3708

marker_gene_table <- marker_gene_table %>% 
  filter(marker_gene_table$gene %in% cell.surface.marker)
saveRDS(marker_gene_table,paste0(outdir,"/hlcma0001/non_stratify_marker_selection_table/old/",split[k],"_SurfaceGenie_marker_gene_table_all.rds"))

seurat <- getMarkers(
  ob_split4,assay = 'RNA',
  organism = 'hg',
  groups = c('annotation.l2'),
  name = 'cerebro_seurat',
  only_pos = F,
  min_pct = 0.5,
  thresh_logFC = 1,
  thresh_p_val = 0.05,
  test = 'wilcox',
  verbose = TRUE
)
marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]] #3708

marker_gene_table <- marker_gene_table %>% 
  filter(marker_gene_table$gene %in% cell.surface.marker)
saveRDS(marker_gene_table,paste0(outdir,"/hlcma0001/level2_marker_selection_table/old/",split[k],"_SurfaceGenie_marker_gene_table_all.rds"))

}  