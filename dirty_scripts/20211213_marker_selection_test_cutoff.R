# Dummy workflow of testing cutoff parameters for marker selection function ####
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

## cell types grouped by master markers' expression created in 01_marker_expression_jessie.R

genes <- c("CD8A","CD8B","CD3D","CD4","CD79A","NKG7","ACTA2")

i <- 5L
print(filename[i])
ob <- readRDS(paste0(rawdir,"/",filename[i],".rds"))
###Compare markers between/among datasets####

split <- as.character(unique(sort(ob@meta.data$disease)))
ob_split <- SplitObject(ob, split.by = "disease")
dir.create(paste0(outdir,"/0_test/",filename[i]))


#for(k in 3:length(split)){ 
k <- 1L
 print(split[k])
  ob_split2 = ob_split[[split[k]]]

  
  
  seurat <- getMarkers(
    ob_split2,assay = 'RNA',
    organism = 'hg',
    groups = c('annotation.l2'),
    name = 'cerebro_seurat',
    only_pos = F,
    genes=genes,
    min_pct = 0.01,
    thresh_logFC = 0.25,
    thresh_p_val = 0.25,
    test = 'wilcox', #DESeq2 no markers found
    verbose = TRUE
  )
  marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]] #3708
  
  marker_gene_table <- marker_gene_table %>% 
    filter(marker_gene_table$gene %in% genes)
  saveRDS(marker_gene_table,paste0(outdir,"/0_test/",filename[i],"/",split[k],"_marker_gene_table_all.rds"))
  
#}  


