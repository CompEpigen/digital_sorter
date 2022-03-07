# Workflow of adding splits AND marker selection of our own lung patient data ####
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(data.table)

# Set working directory ####
wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)
rawdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/patient_hlcma"

outdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/patient_hlcma_addsplits/hlcma0001/3groups_level2"
source("dirty_scripts/getMarkers.R")


cell.surface.marker <-readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/SurfaceGenie/Union.cell.surface.marker.rds")


ob <-readRDS(paste0(rawdir,"/patient_hlcma_tumor_added_addsplits_fixedAnno.rds"))
sampledata <- ob@meta.data
sample_type <- unique(sort(sampledata$disease))
sample_type <- sample_type[c(1,4,3)]

#add annotation.l3####
cellType <- as.character(unique(sort(ob@meta.data$annotation.l2)))
set1 <- c("C1_Tumor","C4_Tumor","C5_Tumor","C15_Tumor")
set2 <- c("Signaling Alveolar Epithelial Type 2","Alveolar Epithelial Type 2" ,
          "Club","Alveolar Epithelial Type 1")
set3 <- c("Mucous" ,"Goblet","Basal","Differentiating Basal","Proliferating Basal",
          "Proximal Basal", 
          "Ciliated" ,"Proximal Ciliated","Ionocyte","Neuroendocrine"  ,"Serous")
meta = ob$annotation.l2
meta[meta %in% set1] = "Tumor"
meta[meta %in% set2] = "Alveolar_Epithelial"
meta[meta %in% set3] = "Airway_Epithelial"
ob <- AddMetaData(object = ob, metadata = meta, col.name = "annotation.l3")


meta = ob$annotation.l3
meta[meta %in% c("Tumor","Alveolar_Epithelial")] = "VS1"
ob <- AddMetaData(object = ob, metadata = meta, col.name = "VS1")

meta = ob$annotation.l3
meta[meta %in% c("Tumor","Airway_Epithelial")] = "VS2"
ob <- AddMetaData(object = ob, metadata = meta, col.name = "VS2")

ob_split <- SplitObject(ob, split.by = "disease")


k <- 3L
for(k in 1:length(sample_type)){ 

print(sample_type[k])
ob_split2 = ob_split[[sample_type[k]]]
ob_split2 <- SplitObject(ob_split2, split.by = "split1")
ob_split2 = ob_split2$`PTPRC-`
ob_split2 <- SplitObject(ob_split2, split.by = "split2")
ob_split2 = ob_split2$`EPCAM+`

ob_split3 <- SplitObject(ob_split2, split.by = "VS1")
ob_split3 = ob_split3[["VS1"]]

ob_split4 <- SplitObject(ob_split2, split.by = "VS2")
ob_split4 = ob_split4[["VS2"]]

seurat <- getMarkers(
  ob_split3,assay = 'RNA',
  organism = 'hg',
  groups = c('annotation.l3'),
  name = 'cerebro_seurat',
  only_pos = F,
  min_pct = 0.5,
  thresh_logFC = 1,
  thresh_p_val = 0.05,
  test = 'wilcox',
  verbose = TRUE
)
marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l3"]] #3708

marker_gene_table <- marker_gene_table %>% 
  filter(marker_gene_table$gene %in% cell.surface.marker)
saveRDS(marker_gene_table,paste0(outdir,"/old/",sample_type[k],"_marker_gene_table_Tumor_Vs_Alveolar_Epithelial.rds"))
if(!is.null(marker_gene_table)){
  # step 0: filter for FDR < 2 (to collect almost all)
  marker_gene_table <- marker_gene_table %>%
    filter(marker_gene_table$p_val_adj<2)
  
    features <- unique(sort(marker_gene_table$gene))
    d <- DotPlot(ob_split3, features = features, group.by = "annotation.l3") + RotatedAxis()
    exp_data <- d[["data"]]
    exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
    exp_data <- exp_data[,c(1,2,5,6)]
    marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l3)
    marker_gene_table <- merge(marker_gene_table, exp_data, by = "mergeID")
    marker_gene_table$backgroud_exp <- marker_gene_table$avg.exp/(2^(marker_gene_table$avg_log2FC))
    
  marker_gene_table <- marker_gene_table %>% 
    arrange(annotation.l3,p_val)
  marker_gene_table$mergeID <- NULL
  
  subset_marker_gene <- data.table(marker_gene_table, key="annotation.l3") 
  # step 0: filter for non-HLA
  subset_marker_gene = subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] 
  # subset_marker_gene <- subset_marker_gene[, head(.SD, 100), by=annotation.l3] 
  print(dim(subset_marker_gene))
}
saveRDS(subset_marker_gene,paste0(outdir,"/add_background/",sample_type[k],"_marker_gene_table_Tumor_Vs_Alveolar_background.rds"))


seurat <- getMarkers(
  ob_split4,assay = 'RNA',
  organism = 'hg',
  groups = c('annotation.l3'),
  name = 'cerebro_seurat',
  only_pos = F,
  min_pct = 0.5,
  thresh_logFC = 1,
  thresh_p_val = 0.05,
  test = 'wilcox',
  verbose = TRUE
)
marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l3"]] #3708

marker_gene_table <- marker_gene_table %>% 
  filter(marker_gene_table$gene %in% cell.surface.marker)
saveRDS(marker_gene_table,paste0(outdir,"/old/",sample_type[k],"_marker_gene_table_Tumor_Vs_Airway_Epithelial.rds"))


if(!is.null(marker_gene_table)){
  # step 0: filter for FDR < 2 (to collect almost all)
  marker_gene_table <- marker_gene_table %>%
    filter(marker_gene_table$p_val_adj<2)
  
  features <- unique(sort(marker_gene_table$gene))
  d <- DotPlot(ob_split4, features = features, group.by = "annotation.l3") + RotatedAxis()
  exp_data <- d[["data"]]
  exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
  exp_data <- exp_data[,c(1,2,5,6)]
  marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l3)
  marker_gene_table <- merge(marker_gene_table, exp_data, by = "mergeID")
  marker_gene_table$backgroud_exp <- marker_gene_table$avg.exp/(2^(marker_gene_table$avg_log2FC))
  
  marker_gene_table <- marker_gene_table %>% 
    arrange(annotation.l3,p_val)
  marker_gene_table$mergeID <- NULL
  
  subset_marker_gene <- data.table(marker_gene_table, key="annotation.l3") 
  # step 0: filter for non-HLA
  subset_marker_gene = subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] 
  # subset_marker_gene <- subset_marker_gene[, head(.SD, 100), by=annotation.l3] 
  print(dim(subset_marker_gene))
}
saveRDS(subset_marker_gene,paste0(outdir,"/add_background/",sample_type[k],"_marker_gene_table_Tumor_Vs_Airway_background.rds"))


}  

