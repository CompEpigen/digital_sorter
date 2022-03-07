# save marker selection table to csv ####
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(data.table)
#install.packages('Seurat')
library(Seurat)
# Set working directory ####
wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)
rawdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/patient_hlcma"

outdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/patient_hlcma_addsplits"


# Read files ####

ob <- readRDS(paste0("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/patient_hlcma/patient_hlcma_tumor_added_addsplits.rds"))

sampledata <- ob@meta.data

sample_type <- unique(sort(sampledata$disease))
ob_split <- SplitObject(ob, split.by = "disease")
## non-stratified cells ####
k <- 2L
for(k in 1:length(sample_type)){ 
  
  print(sample_type[k])
  ob_split2 <- ob_split[[sample_type[k]]]
  if(k==1){
    marker_gene_table <- readRDS(paste0(outdir,"/hlcma0001/non_stratify_marker_selection_table/old/",sample_type[k],"_marker_gene_table_all.rds"))
    
    marker_gene_table$sample_type <- sample_type[k]
    features <- unique(sort(marker_gene_table$gene))
    d <- DotPlot(ob_split2, features = features, group.by = "annotation.l2") + RotatedAxis()
    exp_data <- d[["data"]]
    exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
    exp_data <- exp_data[,c(1,2,5,6)]
    marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
    marker_gene_table <- merge(marker_gene_table, exp_data, by = "mergeID")
    marker_gene_table$backgroud_exp <- marker_gene_table$avg.exp/(2^(marker_gene_table$avg_log2FC))
    #marker_gene_table$mergeID <- NULL
  }else{
    marker_gene_table2 <- readRDS(paste0(outdir,"/hlcma0001/non_stratify_marker_selection_table/old/",sample_type[k],"_marker_gene_table_all.rds"))
    
    marker_gene_table2$sample_type <- sample_type[k]
    features <- unique(sort(marker_gene_table2$gene))
    d <- DotPlot(ob_split2, features = features, group.by = "annotation.l2") + RotatedAxis()
    exp_data <- d[["data"]]
    exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
    exp_data <- exp_data[,c(1,2,5,6)]
    marker_gene_table2$mergeID <- paste(marker_gene_table2$gene, marker_gene_table2$annotation.l2)
    marker_gene_table2 <- merge(marker_gene_table2, exp_data, by = "mergeID")
    marker_gene_table2$backgroud_exp <- marker_gene_table2$avg.exp/(2^(marker_gene_table2$avg_log2FC))
    
    marker_gene_table <- rbind(marker_gene_table,marker_gene_table2)
    
    print(dim(marker_gene_table))
  }
  marker_gene_table$mergeID <- NULL
  write.csv(marker_gene_table,paste0(outdir,"/hlcma0001/non_stratify_marker_selection_table/old/non_stratified_marker_slection_table_expression_added.csv"),row.names = F)
  
}  

## level2 ####
for(k in 1:length(sample_type)){ 
  
  print(sample_type[k])
  ob_split2 <- ob_split[[sample_type[k]]]
  if(k==1){
    marker_gene_table <- readRDS(paste0(outdir,"/hlcma0001/level2_marker_selection_table/old/",sample_type[k],"_marker_gene_table_all.rds"))
    
    marker_gene_table$sample_type <- sample_type[k]
    features <- unique(sort(marker_gene_table$gene))
    d <- DotPlot(ob_split2, features = features, group.by = "annotation.l2") + RotatedAxis()
    exp_data <- d[["data"]]
    exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
    exp_data <- exp_data[,c(1,2,5,6)]
    marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
    marker_gene_table <- merge(marker_gene_table, exp_data, by = "mergeID")
    marker_gene_table$backgroud_exp <- marker_gene_table$avg.exp/(2^(marker_gene_table$avg_log2FC))
    #marker_gene_table$mergeID <- NULL
  }else{
    marker_gene_table2 <- readRDS(paste0(outdir,"/hlcma0001/level2_marker_selection_table/old/",sample_type[k],"_marker_gene_table_all.rds"))
    
    marker_gene_table2$sample_type <- sample_type[k]
    features <- unique(sort(marker_gene_table2$gene))
    d <- DotPlot(ob_split2, features = features, group.by = "annotation.l2") + RotatedAxis()
    exp_data <- d[["data"]]
    exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
    exp_data <- exp_data[,c(1,2,5,6)]
    marker_gene_table2$mergeID <- paste(marker_gene_table2$gene, marker_gene_table2$annotation.l2)
    marker_gene_table2 <- merge(marker_gene_table2, exp_data, by = "mergeID")
    marker_gene_table2$backgroud_exp <- marker_gene_table2$avg.exp/(2^(marker_gene_table2$avg_log2FC))
    
    marker_gene_table <- rbind(marker_gene_table,marker_gene_table2)
    
    print(dim(marker_gene_table))
  }
  marker_gene_table$mergeID <- NULL
  write.csv(marker_gene_table,paste0(outdir,"/hlcma0001/level2_marker_selection_table/old/level2_marker_slection_table_expression_added.csv"),row.names = F)
  
}  


## CD45- ####
for(k in 1:length(sample_type)){ 
  
  print(sample_type[k])
  ob_split2 <- ob_split[[sample_type[k]]]
  if(k==1){
    marker_gene_table <- readRDS(paste0(outdir,"/hlcma0001/CD45-_marker_selection_table/old/",sample_type[k],"_marker_gene_table_all.rds"))
    
    marker_gene_table$sample_type <- sample_type[k]
    features <- unique(sort(marker_gene_table$gene))
    d <- DotPlot(ob_split2, features = features, group.by = "annotation.l2") + RotatedAxis()
    exp_data <- d[["data"]]
    exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
    exp_data <- exp_data[,c(1,2,5,6)]
    marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
    marker_gene_table <- merge(marker_gene_table, exp_data, by = "mergeID")
    marker_gene_table$backgroud_exp <- marker_gene_table$avg.exp/(2^(marker_gene_table$avg_log2FC))
    print(dim(marker_gene_table))
    #marker_gene_table$mergeID <- NULL
  }else{
    marker_gene_table2 <- readRDS(paste0(outdir,"/hlcma0001/CD45-_marker_selection_table/old/",sample_type[k],"_marker_gene_table_all.rds"))
    print(dim(marker_gene_table2))
    marker_gene_table2$sample_type <- sample_type[k]
    features <- unique(sort(marker_gene_table2$gene))
    d <- DotPlot(ob_split2, features = features, group.by = "annotation.l2") + RotatedAxis()
    exp_data <- d[["data"]]
    exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
    exp_data <- exp_data[,c(1,2,5,6)]
    marker_gene_table2$mergeID <- paste(marker_gene_table2$gene, marker_gene_table2$annotation.l2)
    marker_gene_table2 <- merge(marker_gene_table2, exp_data, by = "mergeID")
    marker_gene_table2$backgroud_exp <- marker_gene_table2$avg.exp/(2^(marker_gene_table2$avg_log2FC))
    
    marker_gene_table <- rbind(marker_gene_table,marker_gene_table2)
    
    print(dim(marker_gene_table))
  }
  marker_gene_table$mergeID <- NULL
  write.csv(marker_gene_table,paste0(outdir,"/hlcma0001/CD45-_marker_selection_table/old/CD45-_marker_slection_table_expression_added.csv"),row.names = F)
  
}  
