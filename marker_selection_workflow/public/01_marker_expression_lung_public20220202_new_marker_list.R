#  workflow of marker selection function ####
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(cerebroApp)

# Set working directory ####
wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)
rawdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_each_dataset"
#dir.create(paste0(rawdir,"_addsplits"))
outdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/lung_mapped_cellxgene_fixed_addsplits"
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


#cell.surface.marker <-readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/digital_sorter/R/cell.surface.marker.rds")
cell.surface.marker <-readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/SurfaceGenie/Union.cell.surface.marker.rds")


i <- 1L
for(i in 3:length(filename)){ 
ob <- readRDS(paste0(rawdir,"/",filename[i],".rds"))
###Compare markers between/among datasets####
print(filename[i])
split <- as.character(unique(sort(ob@meta.data$disease)))
ob_split <- SplitObject(ob, split.by = "disease")
#dir.create(paste0(outdir,"/",filename[i]))


for(k in 1:length(split)){ 
#k <- 2L
 print(split[k])
  ob_split2 = ob_split[[split[k]]]
  ob_split3 <- SplitObject(ob_split2, split.by = "split1")
  ob_split3 = ob_split3$`PTPRC-`
  ob_split4 <- SplitObject(ob_split3, split.by = "split2")
  ob_split4 = ob_split4$`EPCAM+`
  
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
  if(!(is.null(marker_gene_table))){
  marker_gene_table <- marker_gene_table %>% 
    filter(marker_gene_table$gene %in% cell.surface.marker) }
  saveRDS(marker_gene_table,paste0(outdir,"/",filename[i],"/level2_marker_selection_table/",split[k],"_SurfaceGenie_marker_gene_table_all.rds"))
  
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
  if(!(is.null(marker_gene_table))){
  marker_gene_table <- marker_gene_table %>% 
    filter(marker_gene_table$gene %in% cell.surface.marker)}
  saveRDS(marker_gene_table,paste0(outdir,"/",filename[i],"/old_marker_selection_table/",split[k],"_SurfaceGenie_marker_gene_table_all.rds"))
}

}


