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

## cell types grouped by master markers' expression created in public/01_marker_expression_jessie.R
normal <- read.csv("procdata/maker_gene_expression_in_normal_lung_add_groups_fixed.csv",row.names = 1)

cell.surface.marker <-readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/SurfaceGenie/Union.cell.surface.marker.rds")


if(F){
  # copy file from May ####
  ob = readRDS("/omics/groups/OE0219/internal/MJMC/P01_NSCLC/P01.2_enriched_cell_components/analysis/azimuth-meta-analysis/human_lung/cellxgene/lung_mapped_cellxgene_tumor_fixedCellAnno.rds")
  table(ob@meta.data$disease)
  #distal_normal general_normal     luadc       tan 
  #  15733          65662          14696        18485
  # distal_normal general_normal          luadc            tan 
  #  15721          65662          14518          18363 
  table(ob@meta.data$donor)
  #   1     2     3  HLD1 HLD20 HLD21 HLD25 HLD26 HLD35  HLD5 
  #9744 28793 27125  5632  8733  4250  8047  7241  7812  7199 
  #    1     2     3  HLD1 HLD20 HLD21 HLD25 HLD26 HLD35  HLD5 
  #9744 28793 27125  5630  8721  4207  8036  7208  7662  7138 
  table(ob@meta.data$dataset_origin)
  #hlcma0001           travaglini_2020 
  # 48914 -> 48602          65662 
  sampledata <- ob@meta.data
  dataset <- unique(sort(sampledata$dataset_origin))
  sampletype <- unique(sort(sampledata$disease))
  
  cell_type <- unique(sort(sampledata$"annotation.l2"))
  new_cell_type <- cell_type[!(cell_type%in% normal$Cell.type)]
  
  #manually add group for new_cell_type (PTPRC- EPCAM+)
  template <- normal[1,]
  C1_Tumor <- c("C1_Tumor","tumor",NA,NA,NA,NA,"PTPRC-","EPCAM+",NA,NA)
  C15_Tumor <- c("C15_Tumor","tumor",NA,NA,NA,NA,"PTPRC-","EPCAM+",NA,NA)
  C4_Tumor <- c("C4_Tumor","tumor",NA,NA,NA,NA,"PTPRC-","EPCAM+",NA,NA)
  C5_Tumor <- c("C5_Tumor","tumor",NA,NA,NA,NA,"PTPRC-","EPCAM+",NA,NA)
  C9_Unknown <- c("C9_Unknown","tumor",NA,NA,NA,NA,"PTPRC-","EPCAM+",NA,NA)
  normal <- rbind(normal,C1_Tumor,C15_Tumor,C4_Tumor,C5_Tumor,C9_Unknown)
  
  # add split ####
  ## CD45 (PTPRC) ####
  features=c("PTPRC")
  
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
  
  set1 <- unique(sort(normal[normal$MME_group =="MME+",]$Cell.type))
  set2 <- unique(sort(normal[normal$MME_group =="MME-",]$Cell.type))
  meta = ob$annotation.l2
  meta[meta %in% set1] = paste0(features, "+")
  meta[meta %in% set2] = paste0(features, "-")
  
  ob <- AddMetaData(object = ob, metadata = meta, col.name = "split4")
  
  
  
  #saveRDS(ob,"/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/patient_hlcma/patient_hlcma_tumor_added_addsplits.rds")
  saveRDS(ob,"/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/patient_hlcma/patient_hlcma_tumor_added_addsplits_fixedAnno.rds")
  
}



ob <-readRDS(paste0(rawdir,"/patient_hlcma_tumor_added_addsplits_fixedAnno.rds"))

###Compare markers between/among datasets####
split <- as.character(unique(sort(ob@meta.data$disease)))
ob_split <- SplitObject(ob, split.by = "disease")


k <- 2L
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
saveRDS(marker_gene_table,paste0(outdir,"/hlcma0001/non_stratify_marker_selection_table/old/",split[k],"_marker_gene_table_all.rds"))

seurat <- getMarkers(
  ob_split3,assay = 'RNA',
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
saveRDS(marker_gene_table,paste0(outdir,"/hlcma0001/CD45-_marker_selection_table/old/",split[k],"_marker_gene_table_all.rds"))

  
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
saveRDS(marker_gene_table,paste0(outdir,"/hlcma0001/level2_marker_selection_table/old/",split[k],"_marker_gene_table_all.rds"))

}  
