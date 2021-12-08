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
rawdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata"
#dir.create(paste0(rawdir,"_addsplits"))
outdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits"

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
filename <- filename[grepl("*fixed_addsplits",filename)]
## cell types grouped by master markers' expression created in 01_marker_expression_jessie.R
cell_group <- read.csv("procdata/maker_gene_expression_in_normal_lung_add_groups.csv")
  
cell.surface.marker <-readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/digital_sorter/R/cell.surface.marker.rds")
ob <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/lung_mapped_cellxgene_fixed_addsplits.rds")

#dir.create(paste0(outdir,"/",filename))

# Check how many cohorts and sample types in the ob ####
## cohorts ####
table(ob@meta.data$dataset_origin)

## sample types ####
table(ob@meta.data$disease)


###Compare markers between/among datasets####
print("Compare markers between/among datasets")
split <- as.character(unique(sort(ob@meta.data$dataset_origin)))
ob_split <- SplitObject(ob, split.by = "dataset_origin")
#k <- 1L
for(k in 1:length(split)){
  dir.create(paste0(outdir,"/",filename,"/test/",split[k]))
  ob_split2 = ob_split[[split[k]]]
  
  seurat <- cerebroApp::getMarkerGenes(
    ob_split2,assay = 'RNA',
    organism = 'hg',
    groups = c('annotation.l2'),
    name = 'cerebro_seurat',
    only_pos = F,
    min_pct = 0,
    thresh_logFC = 0.25,
    thresh_p_val = 0.01,
    test = 'wilcox', #DESeq2 no markers found
    verbose = TRUE
  )
  marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]] #3708
  
  marker_gene_table <- marker_gene_table %>% 
    filter(marker_gene_table$gene %in% cell.surface.marker)
  saveRDS(marker_gene_table, paste0(outdir,"/",filename,"/test/",split[k],"/marker_gene_table_all.rds"))
  
  ob_split2 <- SplitObject(ob_split2, split.by = "split1")
  
  ## Level 1 ####
    ob_split_PTPRC = ob_split2$`PTPRC+`
    
#######Goal: to have CD8A CD8 in the marker_gene_table!!!!!!!!!
    #CD8A n_cell_surface =="FALSE"??!!
    seurat <- cerebroApp::getMarkerGenes(
      ob_split_PTPRC,assay = 'RNA',
      organism = 'hg',
      groups = c('annotation.l2'),
      name = 'cerebro_seurat',
      only_pos = F,
      min_pct = 0,
      thresh_logFC = 0.25,
      thresh_p_val = 0.01,
      test = 'wilcox', #DESeq2 no markers found
      verbose = TRUE
    )
    marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]] #3708
   
    marker_gene_table <- marker_gene_table %>% 
    filter(marker_gene_table$gene %in% cell.surface.marker)
    saveRDS(marker_gene_table, paste0(outdir,"/",filename,"/test/",split[k],"/marker_gene_table_CD45+.rds"))
    
    
}     

