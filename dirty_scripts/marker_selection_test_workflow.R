#!/usr/bin/env Rscript
# Dummy workflow of marker selection for lung sc-RNA data ####
## Example commands
# Rscript marker_selection_test_workflow.R /omics/groups/OE0219/internal/MJMC/P01_NSCLC/P01.2_enriched_cell_components/analysis/azimuth-meta-analysis/human_lung/results /omics/groups/OE0219/internal/MJMC/P01_NSCLC/P01.2_enriched_cell_components/analysis/azimuth-meta-analysis/human_lung/results/by_cohorts /omics/groups/OE0219/internal/MJMC/P01_NSCLC/P01.2_enriched_cell_components/analysis/azimuth-meta-analysis/human_lung/results/by_cohorts_addsplits /omics/groups/OE0219/internal/MJMC/P01_NSCLC/P01.4_digital_sorter/digital_sorter/dirty_scripts
args <- commandArgs(trailingOnly = TRUE)
# args = c("/omics/groups/OE0219/internal/MJMC/P01_NSCLC/P01.2_enriched_cell_components/analysis/azimuth-meta-analysis/human_lung/results", 
#          "/omics/groups/OE0219/internal/MJMC/P01_NSCLC/P01.2_enriched_cell_components/analysis/azimuth-meta-analysis/human_lung/results/by_cohorts", 
#          "/omics/groups/OE0219/internal/MJMC/P01_NSCLC/P01.2_enriched_cell_components/analysis/azimuth-meta-analysis/human_lung/results/by_cohorts_addsplits",
#          "/omics/groups/OE0219/internal/MJMC/P01_NSCLC/P01.4_digital_sorter/digital_sorter/dirty_scripts")


# Auto-check package dependancy
# R version: R4.1.0
checkPkg <- function(pkg){
  return(requireNamespace(pkg, quietly = TRUE))
}
if(!checkPkg("BiocManager")) install.packages("BiocManager")
if(!checkPkg("devtools")) install.packages("devtools")
if(!checkPkg("cerebroApp")) BiocManager::install('romanhaa/cerebroApp')

# Load libraries
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(cerebroApp)

# Set working directory ####
wd = args[1]
setwd(wd)
rawdir = args[2]
#dir.create(paste0(rawdir,"_addsplits"))
outdir = args[3]

# Load functions
if(T){ 
ReadRDSFiles <- function(fileDir, envir = .GlobalEnv) {
  
  # pattern for file searching
  p <- ".rds$"
  rds <- list.files(path = fileDir, pattern = p)
  out <- vapply(rds, FUN = function(.x) {
    nm <- sub(pattern = p, replacement = "", x = .x)
    # read data in
    # out <- readRDS(file = x)
    
    # load in global env
    assign(nm, value = readRDS(file = file.path(fileDir, .x)), envir = envir)
    if (!exists(nm, envir = envir)) return(FALSE)
    TRUE
  }, FUN.VALUE = logical(1L), USE.NAMES = FALSE)
  
  if (!all(out)) warning("Some `.rds` files not loaded.", call. = FALSE)
  
  spc <- paste0(rep('*', times = nchar(fileDir) + 1), collapse = "")
  cat("RDS Files loaded from:", fileDir, spc, rds[out], sep = "\n ", spc)
}

GetfileNames <- function(fileDir, pattern = ".rds"){
  filenames <- list.files(fileDir, pattern= pattern, full.names=TRUE)
  filename <- sub(paste0(".*",fileDir,"/"), "", filenames)
  filename <- sub(pattern, "", filename)  
  return(filename)
}
}

# Read files ####
#ReadRDSFiles(rawdir, envir = .GlobalEnv)
filename <- GetfileNames(rawdir, pattern = ".rds")
filename <- filename[!grepl("*_addsplits",filename)]
## cell types grouped by master markers' expression created in 01_marker_expression_jessie.R
cell_group <- read.csv(file.path(args[4],"maker_gene_expression_in_normal_lung_add_groups.csv"))
  


#i <- 2L
for(i in 1:length(filename)){
  ob = readRDS(paste0(rawdir,"/",filename[i],".rds"))
  cols = c("steelblue","darkred","gold","coral2")
# Categorizing cells with markers ####
  ## CD45 (PTPRC) ####
  
  features=c("PTPRC")
  m1 =VlnPlot(ob, features, split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
  m1
  
  d = m1[[1]]$data
  
  
  set1 <- unique(sort(cell_group[cell_group$PTPRC_group =="PTPRC+",]$Cell.type))
  set2 <- unique(sort(cell_group[cell_group$PTPRC_group =="PTPRC-",]$Cell.type))
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
  
  
  set1 <- unique(sort(cell_group[cell_group$EPCAM_group =="EPCAM+",]$Cell.type))
  set2 <- unique(sort(cell_group[cell_group$EPCAM_group =="EPCAM-",]$Cell.type))
  meta = ob$annotation.l2
  meta[meta %in% set1] = paste0(features, "+")
  meta[meta %in% set2] = paste0(features, "-")
  
  ob <- AddMetaData(object = ob, metadata = meta, col.name = "split2")
  ob2 <- AddMetaData(object = ob2, metadata = meta, col.name = "split2")
  
  
  ## PECAM1####
  ob3 <- SplitObject(ob2, split.by = "split2")
  ob3 = ob3$`EPCAM-`
  
  features=c("PECAM1")
  m1 =VlnPlot(ob3, features, split.by = "disease", group.by = "annotation.l2", 
              cols = cols,sort = TRUE,pt.size = 0, combine = FALSE)
  
  m1
  d = m1[[1]]$data
  
  set1 <- unique(sort(cell_group[cell_group$PECAM1_group =="PECAM1+",]$Cell.type))
  set2 <- unique(sort(cell_group[cell_group$PECAM1_group =="PECAM1-",]$Cell.type))
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
  
  m1
  d = m1[[1]]$data
  
  set1 <- unique(sort(cell_group[cell_group$MME_group =="MME+",]$Cell.type))
  set2 <- unique(sort(cell_group[cell_group$MME_group =="MME-",]$Cell.type))
  meta = ob$annotation.l2
  meta[meta %in% set1] = paste0(features, "+")
  meta[meta %in% set2] = paste0(features, "-")
  
  ob <- AddMetaData(object = ob, metadata = meta, col.name = "split4")
  
  saveRDS(ob,paste0(outdir,"/",filename[i],"_addsplits.rds"))
  
  dir.create(paste0(outdir,"/",filename[i],"_addsplits"))
  
  # Check how many cohorts and sample types in the ob ####
  ## cohorts ####
  table(ob@meta.data$dataset_origin)
  num_cohort <- length(unique(sort(ob@meta.data$dataset_origin)))
  ## sample types ####
  table(ob@meta.data$disease)
  num_sample_type <- length(unique(sort(ob@meta.data$disease)))
  
  if(num_cohort>1){
    ###Compare markers between/among datasets####
    print("Compare markers between/among datasets")
    split <- as.character(unique(sort(ob@meta.data$dataset_origin)))
    ob_split <- SplitObject(ob, split.by = "dataset_origin")
    #k <- 1L
    for(k in 1:length(split)){
      dir.create(paste0(outdir,"/",filename[i],"_addsplits/",split[k]))
      ob_split2 = ob_split[[split[k]]]
      ob_split2 <- SplitObject(ob_split2, split.by = "split1.1")
      ## Level 1 ####
      ob_split_PTPRC = ob_split2$`PTPRC+`
      
      seurat <- cerebroApp::getMarkerGenes(
        ob_split_PTPRC,assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = T,
        min_pct = 0.7,
        thresh_logFC = 0.25,
        thresh_p_val = 0.01,
        test = 'wilcox', #DESeq2 no markers found
        verbose = TRUE
      )
      marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]] #3708
      marker_gene_table <- marker_gene_table %>% 
        filter(marker_gene_table$on_cell_surface =="TRUE")
      
      # step 0: filter for FDR < 0.05
      marker_gene_table <- marker_gene_table %>%
        filter(marker_gene_table$p_val_adj<0.05)%>% 
        arrange(desc(avg_log2FC))
      
      features <- unique(sort(marker_gene_table$gene))
      # plot the initial dot plot
      pdf(paste0(outdir,"/",filename[i],"_addsplits/",split[k],"/Dotplot_level1.pdf"), width=14, height=6)
      d <- DotPlot(ob_split_PTPRC, features = features, group.by = "annotation.l2") + RotatedAxis()
      d
      dev.off()
      # step 1: filter for the average expression > 1
      exp_data <- d[["data"]]
      exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
      exp_data <- exp_data[,c(1,2,5,6)]
      marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
      marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
      
      marker_gene_table2 <- marker_gene_table2 %>%
        filter(marker_gene_table2$avg.exp >1)   
      
      # step 1.5: filter out frequently appearing markers 
      num <- length(table(marker_gene_table2$annotation.l2))/2+1
      repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <=num )
      colnames(repeat_genes) <- "repeat_genes"
      repeat_genes$gene <- row.names(repeat_genes)
      marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
      marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] 
      
      # step 2: arrange according to abs(log2FC)
      marker_gene_table3 <- marker_gene_table3 %>% 
        arrange(annotation.l2,desc(abs(avg_log2FC)))
      marker_gene_table3$mergeID <- NULL
      
      # step 3: select top 10 for each cell type (exclude HLA related)
      library(data.table)
      subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") 
      subset_marker_gene=subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] 
      subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2]
      subset_marker_gene$repeat_genes <- NULL
      selected_marker <- subset_marker_gene
      
      ## Level 2 ####
      ob2 = ob_split2$`PTPRC-`
      ob22 <- SplitObject(ob2, split.by = "split2")
      ob_split_PTPRC0 = ob22$`EPCAM+`
      
      seurat <- cerebroApp::getMarkerGenes(
        ob_split_PTPRC0,assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = T,
        min_pct = 0.7,
        thresh_logFC = 0.25,
        thresh_p_val = 0.01,
        test = 'wilcox', #DESeq2 no markers found
        verbose = TRUE
      )
      marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
      marker_gene_table <- marker_gene_table %>% 
        filter(marker_gene_table$on_cell_surface =="TRUE") 
      
      # step 0: filter for FDR < 0.05
      marker_gene_table <- marker_gene_table %>%
        filter(marker_gene_table$p_val_adj<0.05)%>% 
        arrange(desc(avg_log2FC))
      features <- unique(sort(marker_gene_table$gene))
      
      # plot the initial dot plot
      pdf(paste0(outdir,"/",filename[i],"_addsplits/",split[k],"/Dotplot_level2.pdf"), width=20, height=6)
      d <- DotPlot(ob_split_PTPRC0, features = features, group.by = "annotation.l2") + RotatedAxis()
      d
      dev.off()
      
      # step 1: filter for the average expression > 1
      exp_data <- d[["data"]]
      exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
      exp_data <- exp_data[,c(1,2,5,6)]
      marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
      marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
      
      marker_gene_table2 <- marker_gene_table2 %>%
        filter(marker_gene_table2$avg.exp >1)   
      
      # step 1.5: filter out frequently appearing markers
      table(marker_gene_table2$annotation.l2)
      num <- length(table(marker_gene_table2$annotation.l2))/2+1
      table(marker_gene_table2$gene)
      repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <=num )
      colnames(repeat_genes) <- "repeat_genes"
      repeat_genes$gene <- row.names(repeat_genes)
      marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
      marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] #106
      
      # step 2: arrange according to abs(log2FC)
      marker_gene_table3 <- marker_gene_table3 %>% 
        arrange(annotation.l2,desc(abs(avg_log2FC)))
      marker_gene_table3$mergeID <- NULL
      
      # step 3: select top 10 for each cell type
      library(data.table)
      subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") 
      subset_marker_gene=subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] 
      subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] 
      subset_marker_gene$repeat_genes <- NULL
      selected_marker <- rbind(selected_marker,subset_marker_gene)
      
      
      ## Level 3 ####
      ob3 = ob22$`EPCAM-`
      ob3 <- SplitObject(ob3, split.by = "split3")
      ob_split_EPCAM0 = ob3$`PECAM1+` #EPCAM-
      seurat <- cerebroApp::getMarkerGenes(
        ob_split_EPCAM0,assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = T,
        min_pct = 0.7,
        thresh_logFC = 0.25,
        thresh_p_val = 0.01,
        test = 'wilcox', 
        verbose = TRUE
      )
      marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
      marker_gene_table <- marker_gene_table %>% 
        filter(marker_gene_table$on_cell_surface =="TRUE") #56 #118
      
      # step 0: filter for FDR < 0.05
      marker_gene_table <- marker_gene_table %>%
        filter(marker_gene_table$p_val_adj<0.05)%>% #5 #113
        arrange(desc(avg_log2FC)) 
      
      features <- unique(sort(marker_gene_table$gene))
      # plot the initial dot plot
      pdf(paste0(outdir,"/",filename[i],"_addsplits/",split[k],"/Dotplot_level3.pdf"), width=11, height=6)
      d <- DotPlot(ob_split_EPCAM0, features = features, group.by = "annotation.l2") + RotatedAxis()
      d
      dev.off()
      
      # step 1: filter for the average expression > 1
      exp_data <- d[["data"]]
      exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
      exp_data <- exp_data[,c(1,2,5,6)]
      marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
      marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
      
      marker_gene_table2 <- marker_gene_table2 %>%
        filter(marker_gene_table2$avg.exp >1)   #5 #101
      
      # step 1.5: filter out frequently appearing markers
      table(marker_gene_table2$annotation.l2)
      num <- length(table(marker_gene_table2$annotation.l2))/2+1
      table(marker_gene_table2$gene)
      repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <= num )
      colnames(repeat_genes) <- "repeat_genes"
      repeat_genes$gene <- row.names(repeat_genes)
      marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
      marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] #65 #101
      
      # step 2: arrange according to abs(log2FC)
      marker_gene_table3 <- marker_gene_table3 %>% 
        arrange(annotation.l2,desc(abs(avg_log2FC)))
      
      marker_gene_table3$mergeID <- NULL
      # step 3: select top 10 for each cell type
      library(data.table)
      subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #5
      subset_marker_gene <- subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] #5 #71
      subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] #45 #53
      subset_marker_gene$repeat_genes <- NULL
      selected_marker <- rbind(selected_marker,subset_marker_gene)
      
      ## Level 4 ####
      ob_split_PECAM1 = ob3$`PECAM1-`
      ob4 <- SplitObject(ob_split_PECAM1, split.by = "split4")
      seurat <- cerebroApp::getMarkerGenes(
        ob_split_PECAM1,assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = T,
        min_pct = 0.7,
        thresh_logFC = 0.25,
        thresh_p_val = 0.01,
        test = 'wilcox',
        verbose = TRUE
      )
      marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
      marker_gene_table <- marker_gene_table %>% 
        filter(marker_gene_table$on_cell_surface =="TRUE") #57 #171
      
      # step 0: filter for FDR < 0.05
      marker_gene_table <- marker_gene_table %>%
        filter(marker_gene_table$p_val_adj<0.05)%>% #20 #124
        arrange(desc(avg_log2FC))
      features <- unique(sort(marker_gene_table$gene))
      
      # plot the initial dot plot
      pdf(paste0(outdir,"/",filename[i],"_addsplits/",split[k],"/Dotplot_level4.pdf"), width=18, height=6)
      d <- DotPlot(ob_split_PECAM1, features = features, group.by = "annotation.l2") + RotatedAxis()
      d
      dev.off()
      
      # step 1: filter for the average expression > 1
      exp_data <- d[["data"]]
      exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
      exp_data <- exp_data[,c(1,2,5,6)]
      marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
      marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
      
      marker_gene_table2 <- marker_gene_table2 %>%
        filter(marker_gene_table2$avg.exp >1)   #14 #113
      
      # step 1.5: filter out frequently appearing markers
      num <- length(table(marker_gene_table2$annotation.l2))/2+1
      repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <= num )
      colnames(repeat_genes) <- "repeat_genes"
      repeat_genes$gene <- row.names(repeat_genes)
      marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
      marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] #14 #99
      
      # step 2: arrange according to abs(log2FC)
      marker_gene_table3 <- marker_gene_table3 %>% 
        arrange(annotation.l2,desc(abs(avg_log2FC)))
      marker_gene_table3$mergeID <- NULL
      
      # step 3: select top 10 for each cell type
      library(data.table)
      subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #14 
      subset_marker_gene=subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] #14 #89
      subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] #14 #69
      subset_marker_gene$repeat_genes <- NULL
      selected_marker <- rbind(selected_marker,subset_marker_gene)
      saveRDS(selected_marker,paste0(outdir,"/",filename[i],"_addsplits/",split[k],"/marker_selected_all_levels.rds"))
    print("END--Compare markers between/among datasets--END")
      }
    
    
   
  }
  
  
  if(num_sample_type>1){
    ### Compare markers between/among sample types ####
    print("Compare markers between/among sample types")
    split <- as.character(unique(sort(ob@meta.data$disease)))
    ob_split <- SplitObject(ob, split.by = "disease")
    #k <- 2L
    for(k in 1:length(split)){
      dir.create(paste0(outdir,"/",filename[i],"_addsplits/",split[k]))
      ob_split2 = ob_split[[split[k]]]
      ob_split2 <- SplitObject(ob_split2, split.by = "split1.1")
      ## Level 1 ####
      ob_split_PTPRC = ob_split2$`PTPRC+`
      
      seurat <- cerebroApp::getMarkerGenes(
        ob_split_PTPRC,assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = T,
        min_pct = 0.7,
        thresh_logFC = 0.25,
        thresh_p_val = 0.01,
        test = 'wilcox', #DESeq2 no markers found
        verbose = TRUE
      )
      marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]] 
      marker_gene_table <- marker_gene_table %>% 
        filter(marker_gene_table$on_cell_surface =="TRUE")#279 #266
      
      # step 0: filter for FDR < 0.05
      marker_gene_table <- marker_gene_table %>%
        filter(marker_gene_table$p_val_adj<0.05)%>% #222 #184
        arrange(desc(avg_log2FC))
      
      features <- unique(sort(marker_gene_table$gene))
      # plot the initial dot plot
      pdf(paste0(outdir,"/",filename[i],"_addsplits/",split[k],"/Dotplot_level1.pdf"), width=14, height=6)
      d <- DotPlot(ob_split_PTPRC, features = features, group.by = "annotation.l2") + RotatedAxis()
      d
      dev.off()
      # step 1: filter for the average expression > 1
      exp_data <- d[["data"]]
      exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
      exp_data <- exp_data[,c(1,2,5,6)]
      marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
      marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
      
      marker_gene_table2 <- marker_gene_table2 %>%
        filter(marker_gene_table2$avg.exp >1)   #221 #171
      
      # step 1.5: filter out frequently appearing markers 
      num <- length(table(marker_gene_table2$annotation.l2))/2+1
      repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <=num )
      colnames(repeat_genes) <- "repeat_genes"
      repeat_genes$gene <- row.names(repeat_genes)
      marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
      marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] #210 #160
      
      # step 2: arrange according to abs(log2FC)
      marker_gene_table3 <- marker_gene_table3 %>% 
        arrange(annotation.l2,desc(abs(avg_log2FC)))
      marker_gene_table3$mergeID <- NULL
      
      # step 3: select top 10 for each cell type (exclude HLA related)
      library(data.table)
      subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #210 #160
      subset_marker_gene=subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] #135 #111
      subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] #111 #102
      subset_marker_gene$repeat_genes <- NULL
      selected_marker <- subset_marker_gene
      
      ## Level 2 ####
      ob2 = ob_split2$`PTPRC-`
      ob22 <- SplitObject(ob2, split.by = "split2")
      ob_split_PTPRC0 = ob22$`EPCAM+`
      
      seurat <- cerebroApp::getMarkerGenes(
        ob_split_PTPRC0,assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = T,
        min_pct = 0.7,
        thresh_logFC = 0.25,
        thresh_p_val = 0.01,
        test = 'wilcox', #DESeq2 no markers found
        verbose = TRUE
      )
      marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
      marker_gene_table <- marker_gene_table %>% 
        filter(marker_gene_table$on_cell_surface =="TRUE") #262 #172
      
      # step 0: filter for FDR < 0.05
      marker_gene_table <- marker_gene_table %>%
        filter(marker_gene_table$p_val_adj<0.05)%>% #159 #115
        arrange(desc(avg_log2FC))
      features <- unique(sort(marker_gene_table$gene))
      
      # plot the initial dot plot
      pdf(paste0(outdir,"/",filename[i],"_addsplits/",split[k],"/Dotplot_level2.pdf"), width=20, height=6)
      d <- DotPlot(ob_split_PTPRC0, features = features, group.by = "annotation.l2") + RotatedAxis()
      d
      dev.off()
      
      # step 1: filter for the average expression > 1
      exp_data <- d[["data"]]
      exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
      exp_data <- exp_data[,c(1,2,5,6)]
      marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
      marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
      
      marker_gene_table2 <- marker_gene_table2 %>%
        filter(marker_gene_table2$avg.exp >1)   #140 #106
      
      # step 1.5: filter out frequently appearing markers
      table(marker_gene_table2$annotation.l2)
      num <- length(table(marker_gene_table2$annotation.l2))/2+1
      table(marker_gene_table2$gene)
      repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <=num )
      colnames(repeat_genes) <- "repeat_genes"
      repeat_genes$gene <- row.names(repeat_genes)
      marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
      marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] #106
      
      # step 2: arrange according to abs(log2FC)
      marker_gene_table3 <- marker_gene_table3 %>% 
        arrange(annotation.l2,desc(abs(avg_log2FC)))
      marker_gene_table3$mergeID <- NULL
      
      # step 3: select top 10 for each cell type
      library(data.table)
      subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #140 #106
      subset_marker_gene=subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] #118 #87
      subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] #82 #59
      subset_marker_gene$repeat_genes <- NULL
      selected_marker <- rbind(selected_marker,subset_marker_gene)
      
      
      ## Level 3 ####
      ob3 = ob22$`EPCAM-`
      ob3 <- SplitObject(ob3, split.by = "split3")
      ob_split_EPCAM0 = ob3$`PECAM1+` #EPCAM-
      seurat <- cerebroApp::getMarkerGenes(
        ob_split_EPCAM0,assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = T,
        min_pct = 0.7,
        thresh_logFC = 0.25,
        thresh_p_val = 0.01,
        test = 'wilcox', 
        verbose = TRUE
      )
      marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
      marker_gene_table <- marker_gene_table %>% 
        filter(marker_gene_table$on_cell_surface =="TRUE") #22 #56
      
      # step 0: filter for FDR < 0.1 !!!!!!!!
      marker_gene_table <- marker_gene_table %>%
        filter(marker_gene_table$p_val_adj<0.1)%>% #1 #if 0.05 -> 0 left #5 (also when 0.05)
        arrange(desc(avg_log2FC)) 
      
      features <- unique(sort(marker_gene_table$gene))
      # plot the initial dot plot
      pdf(paste0(outdir,"/",filename[i],"_addsplits/",split[k],"/Dotplot_level3.pdf"), width=11, height=6)
      d <- DotPlot(ob_split_EPCAM0, features = features, group.by = "annotation.l2") + RotatedAxis()
      d
      dev.off()
      
      # step 1: filter for the average expression > 1
      exp_data <- d[["data"]]
      exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
      exp_data <- exp_data[,c(1,2,5,6)]
      marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
      marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
      
      marker_gene_table2 <- marker_gene_table2 %>%
        filter(marker_gene_table2$avg.exp >1)   #1 #5
      
      # step 1.5: filter out frequently appearing markers
      table(marker_gene_table2$annotation.l2)
      num <- length(table(marker_gene_table2$annotation.l2))/2+1
      table(marker_gene_table2$gene)
      repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <= num )
      colnames(repeat_genes) <- "repeat_genes"
      repeat_genes$gene <- row.names(repeat_genes)
      marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
      marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] #1 #5
      
      # step 2: arrange according to abs(log2FC)
      marker_gene_table3 <- marker_gene_table3 %>% 
        arrange(annotation.l2,desc(abs(avg_log2FC)))
      
      marker_gene_table3$mergeID <- NULL
      # step 3: select top 10 for each cell type
      library(data.table)
      subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #1
      subset_marker_gene <- subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] #1 #5
      subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] #1 #5
      subset_marker_gene$repeat_genes <- NULL
      selected_marker <- rbind(selected_marker,subset_marker_gene)
      
      ## Level 4 ####
      ob_split_PECAM1 = ob3$`PECAM1-`
      ob4 <- SplitObject(ob_split_PECAM1, split.by = "split4")
      seurat <- cerebroApp::getMarkerGenes(
        ob_split_PECAM1,assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = T,
        min_pct = 0.7,
        thresh_logFC = 0.25,
        thresh_p_val = 0.01,
        test = 'wilcox',
        verbose = TRUE
      )
      marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
      marker_gene_table <- marker_gene_table %>% 
        filter(marker_gene_table$on_cell_surface =="TRUE") #63 #57
      
      # step 0: filter for FDR < 0.05
      marker_gene_table <- marker_gene_table %>%
        filter(marker_gene_table$p_val_adj<0.05)%>% #42 #20
        arrange(desc(avg_log2FC))
      features <- unique(sort(marker_gene_table$gene))
      
      # plot the initial dot plot
      pdf(paste0(outdir,"/",filename[i],"_addsplits/",split[k],"/Dotplot_level4.pdf"), width=18, height=6)
      d <- DotPlot(ob_split_PECAM1, features = features, group.by = "annotation.l2") + RotatedAxis()
      d
      dev.off()
      
      # step 1: filter for the average expression > 1
      exp_data <- d[["data"]]
      exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
      exp_data <- exp_data[,c(1,2,5,6)]
      marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
      marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
      
      marker_gene_table2 <- marker_gene_table2 %>%
        filter(marker_gene_table2$avg.exp >1)   #28 #14
      
      # step 1.5: filter out frequently appearing markers
      num <- length(table(marker_gene_table2$annotation.l2))/2+1
      repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <= num )
      colnames(repeat_genes) <- "repeat_genes"
      repeat_genes$gene <- row.names(repeat_genes)
      marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
      marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] #28 #14
      
      # step 2: arrange according to abs(log2FC)
      marker_gene_table3 <- marker_gene_table3 %>% 
        arrange(annotation.l2,desc(abs(avg_log2FC)))
      marker_gene_table3$mergeID <- NULL
      
      # step 3: select top 10 for each cell type
      library(data.table)
      subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #28 
      subset_marker_gene=subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] #27 #14
      subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] #19 #14
      subset_marker_gene$repeat_genes <- NULL
      selected_marker <- rbind(selected_marker,subset_marker_gene)
      saveRDS(selected_marker,paste0(outdir,"/",filename[i],"_addsplits/",split[k],"/marker_selected_all_levels.rds"))
      print("END--Compare markers between/among datasets--END")
  }
}


 } 


