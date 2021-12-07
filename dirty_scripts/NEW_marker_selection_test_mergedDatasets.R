# Dummy workflow of marker selection for lung sc-RNA data ####
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
if(T){ 

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
filename <- filename[grepl("*gene_addsplits",filename)]

cell.surface.marker <-readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/digital_sorter/R/cell.surface.marker.rds")
ob <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/lung_mapped_cellxgene_addsplits.rds")

# Check how many cohorts and sample types in the ob ####
## cohorts ####
table(ob@meta.data$dataset_origin)
num_cohort <- length(unique(sort(ob@meta.data$dataset_origin)))
## sample types ####
table(ob@meta.data$disease)
num_sample_type <- length(unique(sort(ob@meta.data$disease)))
  
#if(num_cohort>1){
if(F){
  ###Compare markers between/among datasets####
  print("Compare markers between/among datasets")
  split <- as.character(unique(sort(ob@meta.data$dataset_origin)))
  ob_split <- SplitObject(ob, split.by = "dataset_origin")
  k <- 1L
  for(k in 1:length(split)){
    dir.create(paste0(outdir,"/",filename,"/by_datasets/",split[k]))
    ob_split2 = ob_split[[split[k]]]
    ob_split2 <- SplitObject(ob_split2, split.by = "split1")
    
    ## Level 1 ####
    ob_split_PTPRC = ob_split2$`PTPRC+`
      
    seurat <- cerebroApp::getMarkerGenes(
        ob_split_PTPRC,assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = F,
        min_pct = 0.5,
        thresh_logFC = 0.2,
        thresh_p_val = 0.1,
        test = 'wilcox', #DESeq2 no markers found
        verbose = TRUE
    )
    marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]] #3708
    if(length(rownames(marker_gene_table))>0){
      if(length(rownames(marker_gene_table))>0){
      marker_gene_table <- marker_gene_table %>% 
        filter(marker_gene_table$gene %in% cell.surface.marker)
      
      # step 0: filter for FDR < 0.1
      marker_gene_table <- marker_gene_table %>%
        filter(marker_gene_table$p_val_adj<0.1)%>% 
        arrange(desc(avg_log2FC))
      }
      features <- unique(sort(marker_gene_table$gene))
      
      # plot the initial dot plot
      if(length(features)>0){
      pdf(paste0(outdir,"/",filename,"/by_datasets/",split[k],"/Dotplot_level1.pdf"), width=14, height=6)
      DotPlot(ob_split_PTPRC, features = features, group.by = "annotation.l2") + RotatedAxis()
      dev.off()
      
      # step 1: filter for the average expression > 1
      d <- DotPlot(ob_split_PTPRC, features = features, group.by = "annotation.l2") + RotatedAxis()
      exp_data <- d[["data"]]
      exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
      exp_data <- exp_data[,c(1,2,5,6)]
      
      }
      
      marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
      marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
      
      if(length(rownames(marker_gene_table2))>0){
      marker_gene_table2 <- marker_gene_table2 %>%
        filter(marker_gene_table2$avg.exp >1)   
      }
      
      # step 1.5: filter out frequently appearing markers 
      num <- length(table(marker_gene_table2$annotation.l2))/2+1
      repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <=num )
      colnames(repeat_genes) <- "repeat_genes"
      repeat_genes$gene <- row.names(repeat_genes)
      marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
      #marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] 
      
      # step 2: arrange according to abs(log2FC)
      marker_gene_table3 <- marker_gene_table3 %>% 
        arrange(annotation.l2,desc(abs(avg_log2FC)))
      marker_gene_table3$mergeID <- NULL
      
      # step 3: select top 10 for each cell type (exclude HLA related)
      library(data.table)
      subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") 
      subset_marker_gene=subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] 
      subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2]
      
      selected_marker <- subset_marker_gene
    }
    
    ## Level 2 ####
    ob2 = ob_split2$`PTPRC-`
    ob22 <- SplitObject(ob2, split.by = "split2")
    ob_split_PTPRC0 = ob22$`EPCAM+`
      
    seurat <- cerebroApp::getMarkerGenes(
        ob_split_PTPRC0,assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = F,
        min_pct = 0.5,
        thresh_logFC = 0.2,
        thresh_p_val = 0.1,
        test = 'wilcox', #DESeq2 no markers found
        verbose = TRUE
      )
    marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
    if(length(rownames(marker_gene_table))>0){
      if(length(rownames(marker_gene_table))>0){
      marker_gene_table <- marker_gene_table %>% 
        filter(marker_gene_table$gene %in% cell.surface.marker) 
      
      # step 0: filter for FDR < 0.1
      marker_gene_table <- marker_gene_table %>%
        filter(marker_gene_table$p_val_adj<0.1)%>% 
        arrange(desc(avg_log2FC))
      }
      features <- unique(sort(marker_gene_table$gene))
      
      # plot the initial dot plot
      if(length(features)>0){
        pdf(paste0(outdir,"/",filename,"/by_datasets/",split[k],"/Dotplot_level2.pdf"), width=20, height=6)
        DotPlot(ob_split_PTPRC0, features = features, group.by = "annotation.l2") + RotatedAxis()
        dev.off()
        
        
      # step 1: filter for the average expression > 1
        d <- DotPlot(ob_split_PTPRC0, features = features, group.by = "annotation.l2") + RotatedAxis()
     
      exp_data <- d[["data"]]
      exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
      exp_data <- exp_data[,c(1,2,5,6)]
      
      }
      
      marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
      marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
      if(length(rownames(marker_gene_table2))>0){
      marker_gene_table2 <- marker_gene_table2 %>%
        filter(marker_gene_table2$avg.exp >1)   
      }
      
      # step 1.5: filter out frequently appearing markers
      table(marker_gene_table2$annotation.l2)
      num <- length(table(marker_gene_table2$annotation.l2))/2+1
      table(marker_gene_table2$gene)
      repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <=num )
      colnames(repeat_genes) <- "repeat_genes"
      repeat_genes$gene <- row.names(repeat_genes)
      marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
      #marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] #106
      
      # step 2: arrange according to abs(log2FC)
      marker_gene_table3 <- marker_gene_table3 %>% 
        arrange(annotation.l2,desc(abs(avg_log2FC)))
      marker_gene_table3$mergeID <- NULL
      
      # step 3: select top 10 for each cell type
      library(data.table)
      subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") 
      subset_marker_gene=subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] 
      subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] 
      if(length(rownames(subset_marker_gene))>0){
      selected_marker <- rbind(selected_marker,subset_marker_gene)
      }
    }
    
    ## Level 3 ####
    ob3 = ob22$`EPCAM-`
    ob3 <- SplitObject(ob3, split.by = "split3")
    ob_split_EPCAM0 = ob3$`PECAM1+` #EPCAM-
    seurat <- cerebroApp::getMarkerGenes(
        ob_split_EPCAM0,assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = F,
        min_pct = 0.5,
        thresh_logFC = 0.2,
        thresh_p_val = 0.1,
        test = 'wilcox', 
        verbose = TRUE
      )
    marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
    if(length(rownames(marker_gene_table))>0){
      if(length(rownames(marker_gene_table))>0){
      marker_gene_table <- marker_gene_table %>% 
        filter(marker_gene_table$gene %in% cell.surface.marker) #56 #118
      
      # step 0: filter for FDR < 0.1
      marker_gene_table <- marker_gene_table %>%
        filter(marker_gene_table$p_val_adj<0.1)%>% #5 #113
        arrange(desc(avg_log2FC)) 
      }
      features <- unique(sort(marker_gene_table$gene))
      
      # plot the initial dot plot
      if(length(features)>0){
      pdf(paste0(outdir,"/",filename,"/by_datasets/",split[k],"/Dotplot_level3.pdf"), width=11, height=6)
      DotPlot(ob_split_EPCAM0, features = features, group.by = "annotation.l2") + RotatedAxis()
      dev.off()
      
      
      # step 1: filter for the average expression > 1
      d <- DotPlot(ob_split_EPCAM0, features = features, group.by = "annotation.l2") + RotatedAxis()
      
      exp_data <- d[["data"]]
      exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
      exp_data <- exp_data[,c(1,2,5,6)]
      
      }
      
      
      marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
      marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
      if(length(rownames(marker_gene_table2))>0){
      marker_gene_table2 <- marker_gene_table2 %>%
        filter(marker_gene_table2$avg.exp >1)   #5 #101
      }
      
      # step 1.5: filter out frequently appearing markers
      table(marker_gene_table2$annotation.l2)
      num <- length(table(marker_gene_table2$annotation.l2))/2+1
      table(marker_gene_table2$gene)
      repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <= num )
      colnames(repeat_genes) <- "repeat_genes"
      repeat_genes$gene <- row.names(repeat_genes)
      marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
      #marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] #65 #101
      
      # step 2: arrange according to abs(log2FC)
      marker_gene_table3 <- marker_gene_table3 %>% 
        arrange(annotation.l2,desc(abs(avg_log2FC)))
      
      marker_gene_table3$mergeID <- NULL
      # step 3: select top 10 for each cell type
      library(data.table)
      subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #5
      subset_marker_gene <- subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] #5 #71
      subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] #45 #53
      if(length(rownames(subset_marker_gene))>0){
      selected_marker <- rbind(selected_marker,subset_marker_gene)
      }
    }
    ## Level 4 ####
    ob_split_PECAM1 = ob3$`PECAM1-`
   
    seurat <- cerebroApp::getMarkerGenes(
        ob_split_PECAM1,assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = F,
        min_pct = 0.5,
        thresh_logFC = 0.2,
        thresh_p_val = 0.1,
        test = 'wilcox',
        verbose = TRUE
    )
    marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
    
    if(length(rownames(marker_gene_table))>0){
      if(length(rownames(marker_gene_table))>0){
      marker_gene_table <- marker_gene_table %>% 
        filter(marker_gene_table$gene %in% cell.surface.marker) #57 #171
      
      # step 0: filter for FDR < 0.1
      marker_gene_table <- marker_gene_table %>%
        filter(marker_gene_table$p_val_adj<0.1)%>% #20 #124
        arrange(desc(avg_log2FC))
      }
      features <- unique(sort(marker_gene_table$gene))
      
      # plot the initial dot plot
      if(length(features)>0){
      pdf(paste0(outdir,"/",filename,"/by_datasets/",split[k],"/Dotplot_level4.pdf"), width=18, height=6)
      DotPlot(ob_split_PECAM1, features = features, group.by = "annotation.l2") + RotatedAxis()
      dev.off()
      
      
      # step 1: filter for the average expression > 1
      d <- DotPlot(ob_split_PECAM1, features = features, group.by = "annotation.l2") + RotatedAxis()
      exp_data <- d[["data"]]
      exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
      exp_data <- exp_data[,c(1,2,5,6)]
      }
      
      marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
      marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
      if(length(rownames(marker_gene_table2))>0){
      marker_gene_table2 <- marker_gene_table2 %>%
        filter(marker_gene_table2$avg.exp >1)   #14 #113
      }
      # step 1.5: filter out frequently appearing markers
      num <- length(table(marker_gene_table2$annotation.l2))/2+1
      repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <= num )
      colnames(repeat_genes) <- "repeat_genes"
      repeat_genes$gene <- row.names(repeat_genes)
      marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
     # marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] #14 #99
      
      # step 2: arrange according to abs(log2FC)
      marker_gene_table3 <- marker_gene_table3 %>% 
        arrange(annotation.l2,desc(abs(avg_log2FC)))
      marker_gene_table3$mergeID <- NULL
      
      # step 3: select top 10 for each cell type
      library(data.table)
      subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #14 
      subset_marker_gene=subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] #14 #89
      subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] #14 #69
      if(length(rownames(subset_marker_gene))>0){
      selected_marker <- rbind(selected_marker,subset_marker_gene)
      }
    saveRDS(selected_marker,paste0(outdir,"/",filename,"/by_datasets/",split[k],"/marker_selected_all_levels.rds"))
    print(paste0("END-- ",k," dataset --END"))
    }
  }
    print("END--Compare markers between/among datasets--END")
   
}
  
  #}
if(num_sample_type>1){
  ### Compare markers between/among sample types ####
  print("Compare markers between/among sample types")
  split <- as.character(unique(sort(ob@meta.data$disease)))
  ob_split <- SplitObject(ob, split.by = "disease")
  k <- 3L
  for(k in 1:length(split)){
    dir.create(paste0(outdir,"/",filename,"/by_diseases/",split[k]))
    ob_split2 = ob_split[[split[k]]]
    ob_split2 <- SplitObject(ob_split2, split.by = "split1")
      
    ## Level 1 ####
    ob_split_PTPRC = ob_split2$`PTPRC+`
      
      if(is.null(x = ob_split_PTPRC)){ 
        marker_gene_table <- NULL
        print("skip level 1")
      }else{ 
      seurat <- cerebroApp::getMarkerGenes(
        ob_split_PTPRC,assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = F,
        min_pct = 0.5,
        thresh_logFC = 0.2,
        thresh_p_val = 0.1,
        test = 'wilcox', #DESeq2 no markers found
        verbose = TRUE
      )
      
      marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]] 
      }
      
      if(length(rownames(marker_gene_table))>0){
        if(length(rownames(marker_gene_table))>0){
        marker_gene_table <- marker_gene_table %>% 
          filter(marker_gene_table$gene %in% cell.surface.marker)#279 #266
      
      
      # step 0: filter for FDR < 0.1
      marker_gene_table <- marker_gene_table %>%
        filter(marker_gene_table$p_val_adj<0.1)%>% #222 #184
        arrange(desc(avg_log2FC))
      }
      
      features <- unique(sort(marker_gene_table$gene))
      
      # plot the initial dot plot
      if(length(features)>0){
      pdf(paste0(outdir,"/",filename,"/by_diseases/",split[k],"/Dotplot_level1.pdf"), width=20, height=6)
      DotPlot(ob_split_PTPRC, features = features, group.by = "annotation.l2") + RotatedAxis()
      dev.off()
     
      # step 1: filter for the average expression > 1
      d <- DotPlot(ob_split_PTPRC, features = features, group.by = "annotation.l2") + RotatedAxis()
      exp_data <- d[["data"]]
      exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
      exp_data <- exp_data[,c(1,2,5,6)]
      }
      
      marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
      marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
      
      if(length(rownames(marker_gene_table2))>0){
        marker_gene_table2 <- marker_gene_table2 %>%
        filter(marker_gene_table2$avg.exp >1)   #221 #171
      }
      
      # step 1.5: filter out frequently appearing markers 
      num <- length(table(marker_gene_table2$annotation.l2))/2+1
      repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <=num )
      colnames(repeat_genes) <- "repeat_genes"
      repeat_genes$gene <- row.names(repeat_genes)
      marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
      
      #marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] #210 #160
      
      # step 2: arrange according to abs(log2FC)
      marker_gene_table3 <- marker_gene_table3 %>% 
        arrange(annotation.l2,desc(abs(avg_log2FC)))
      marker_gene_table3$mergeID <- NULL
      
      # step 3: select top 10 for each cell type (exclude HLA related)
      library(data.table)
      subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #210 #160
      subset_marker_gene=subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] #135 #111
      subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] #111 #102
      #subset_marker_gene$repeat_genes <- NULL
      selected_marker <- subset_marker_gene
      }
      
    ## Level 2 ####
    ob2 = ob_split2$`PTPRC-`
    ob22 <- SplitObject(ob2, split.by = "split2")
    ob_split_PTPRC0 = ob22$`EPCAM+`
      
    if(is.null(x = ob_split_PTPRC0)){ 
        marker_gene_table <- NULL
        print("skip level 2")
      }else{ 
      seurat <- cerebroApp::getMarkerGenes(
        ob_split_PTPRC0,assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = F,
        min_pct = 0.5,
        thresh_logFC = 0.2,
        thresh_p_val = 0.1,
        test = 'wilcox', #DESeq2 no markers found
        verbose = TRUE
      )
     
      marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
     
      
    if(length(rownames(marker_gene_table))>0){
      if(length(rownames(marker_gene_table))>0){
        marker_gene_table <- marker_gene_table %>% 
        filter(marker_gene_table$gene %in% cell.surface.marker) #262 #172
      
      # step 0: filter for FDR < 0.1
      marker_gene_table <- marker_gene_table %>%
        filter(marker_gene_table$p_val_adj<0.1)%>% #159 #115
        arrange(desc(avg_log2FC))
      }
      features <- unique(sort(marker_gene_table$gene))
      
      # plot the initial dot plot
      if(length(features)>0){
      pdf(paste0(outdir,"/",filename,"/by_diseases/",split[k],"/Dotplot_level2.pdf"), width=20, height=6)
      DotPlot(ob_split_PTPRC0, features = features, group.by = "annotation.l2") + RotatedAxis()
      dev.off()
      
      # step 1: filter for the average expression > 1
      d <- DotPlot(ob_split_PTPRC0, features = features, group.by = "annotation.l2") + RotatedAxis()
      
      exp_data <- d[["data"]]
      exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
      exp_data <- exp_data[,c(1,2,5,6)]
      }
      
      marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
      marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
      if(length(rownames(marker_gene_table2))>0){
      marker_gene_table2 <- marker_gene_table2 %>%
        filter(marker_gene_table2$avg.exp >1)   #140 #106
      }
      # step 1.5: filter out frequently appearing markers
      table(marker_gene_table2$annotation.l2)
      num <- length(table(marker_gene_table2$annotation.l2))/2+1
      table(marker_gene_table2$gene)
      repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <=num )
      colnames(repeat_genes) <- "repeat_genes"
      repeat_genes$gene <- row.names(repeat_genes)
      marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
      #marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] #106
      
      # step 2: arrange according to abs(log2FC)
      marker_gene_table3 <- marker_gene_table3 %>% 
        arrange(annotation.l2,desc(abs(avg_log2FC)))
      marker_gene_table3$mergeID <- NULL
      
      # step 3: select top 10 for each cell type
      library(data.table)
      subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #140 #106
      subset_marker_gene=subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] #118 #87
      subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] #82 #59
      
      if(length(rownames(subset_marker_gene))>0){
      selected_marker <- rbind(selected_marker,subset_marker_gene)
      }
    }
}
      
    ## Level 3 ####
    ob3 = ob22$`EPCAM-`
    ob3 <- SplitObject(ob3, split.by = "split3")
    ob_split_EPCAM0 = ob3$`PECAM1+` #EPCAM-
      
    if(is.null(x = ob_split_EPCAM0)){ 
        print("skip level 3")
      }else{
        seurat <- cerebroApp::getMarkerGenes(
          ob_split_EPCAM0,assay = 'RNA',
          organism = 'hg',
          groups = c('annotation.l2'),
          name = 'cerebro_seurat',
          only_pos = F,
          min_pct = 0.5,
          thresh_logFC = 0.2,
          thresh_p_val = 0.1,
          test = 'wilcox', 
          verbose = TRUE
        )
      marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
      
        if(length(rownames(marker_gene_table))>0){
          if(length(rownames(marker_gene_table))>0){
            marker_gene_table <- marker_gene_table %>% 
              filter(marker_gene_table$gene %in% cell.surface.marker) #22 #56
          }
          # step 0: filter for FDR < 0.1 !!!!!!!!
          marker_gene_table <- marker_gene_table %>%
            filter(marker_gene_table$p_val_adj<0.1)%>% 
            arrange(desc(avg_log2FC)) 
          
          features <- unique(sort(marker_gene_table$gene))
          # plot the initial dot plot
          if(length(features)>0){
            pdf(paste0(outdir,"/",filename,"/by_diseases/",split[k],"/Dotplot_level3.pdf"), width=20, height=6)
            DotPlot(ob_split_EPCAM0, features = features, group.by = "annotation.l2") + RotatedAxis()
        
            dev.off()
          
          # step 1: filter for the average expression > 1
          d <- DotPlot(ob_split_EPCAM0, features = features, group.by = "annotation.l2") + RotatedAxis()
          
          exp_data <- d[["data"]]
          exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
          exp_data <- exp_data[,c(1,2,5,6)]
          }
          marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
          marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
          if(length(rownames(marker_gene_table2))>0){
            marker_gene_table2 <- marker_gene_table2 %>%
              filter(marker_gene_table2$avg.exp >1)   #1 #5
          }
          # step 1.5: filter out frequently appearing markers
          table(marker_gene_table2$annotation.l2)
          num <- length(table(marker_gene_table2$annotation.l2))/2+1
          table(marker_gene_table2$gene)
          repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <= num )
          colnames(repeat_genes) <- "repeat_genes"
          repeat_genes$gene <- row.names(repeat_genes)
          marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
          #marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] #1 #5
          
          # step 2: arrange according to abs(log2FC)
          marker_gene_table3 <- marker_gene_table3 %>% 
            arrange(annotation.l2,desc(abs(avg_log2FC)))
          
          marker_gene_table3$mergeID <- NULL
          # step 3: select top 10 for each cell type
          library(data.table)
          subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #1
          subset_marker_gene <- subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] #1 #5
          subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] #1 #5
          if(length(rownames(subset_marker_gene))>0){
            selected_marker <- rbind(selected_marker,subset_marker_gene)
          }
        }
      }
      
      ## Level 4 ####
      ob_split_PECAM1 = ob3$`PECAM1-`
      #ob4 <- SplitObject(ob_split_PECAM1, split.by = "split4")
      if(is.null(x = ob_split_PECAM1) ){ 
        print("skip level 4")
      }else{
      seurat <- cerebroApp::getMarkerGenes(
        ob_split_PECAM1,assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = F,
        min_pct = 0.5,
        thresh_logFC = 0.2,
        thresh_p_val = 0.1,
        test = 'wilcox',
        verbose = TRUE
      )
      marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
      
    if(length(rownames(marker_gene_table))>0){
      if(length(rownames(marker_gene_table))>0){
      marker_gene_table <- marker_gene_table %>% 
        filter(marker_gene_table$gene %in% cell.surface.marker) #63 #57
      }
      # step 0: filter for FDR < 0.1
      marker_gene_table <- marker_gene_table %>%
        filter(marker_gene_table$p_val_adj<0.1)%>% #42 #20
        arrange(desc(avg_log2FC))
      features <- unique(sort(marker_gene_table$gene))
      
      # plot the initial dot plot
      if(length(features)>0){
      pdf(paste0(outdir,"/",filename,"/by_diseases/",split[k],"/Dotplot_level4.pdf"), width=18, height=6)
      DotPlot(ob_split_PECAM1, features = features, group.by = "annotation.l2") + RotatedAxis()
      dev.off()
      
      # step 1: filter for the average expression > 1
      d <- DotPlot(ob_split_PECAM1, features = features, group.by = "annotation.l2") + RotatedAxis()
      
      exp_data <- d[["data"]]
      exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
      exp_data <- exp_data[,c(1,2,5,6)]
      }
      
      marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
      marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")
      if(length(rownames(marker_gene_table2))>0){
      marker_gene_table2 <- marker_gene_table2 %>%
        filter(marker_gene_table2$avg.exp >1)   #28 #14
      }
      # step 1.5: filter out frequently appearing markers
      num <- length(table(marker_gene_table2$annotation.l2))/2+1
      repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <= num )
      colnames(repeat_genes) <- "repeat_genes"
      repeat_genes$gene <- row.names(repeat_genes)
      marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
      #marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat_genes,] #28 #14
      
      # step 2: arrange according to abs(log2FC)
      marker_gene_table3 <- marker_gene_table3 %>% 
        arrange(annotation.l2,desc(abs(avg_log2FC)))
      marker_gene_table3$mergeID <- NULL
      
      # step 3: select top 10 for each cell type
      library(data.table)
      subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #28 
      subset_marker_gene=subset_marker_gene[!grepl("HLA",subset_marker_gene$gene),] #27 #14
      subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] #19 #14
      if(length(rownames(subset_marker_gene))>0){
      selected_marker <- rbind(selected_marker,subset_marker_gene)
      }
      
      
    }
}
      
      ## Final Save File ####
      saveRDS(selected_marker,paste0(outdir,"/",filename,"/by_diseases/",split[k],"/marker_selected_all_levels.rds"))
      print(paste0("END-- ",k," dataset --END"))
    }
    print("END--Compare markers between/among datasets--END")
}


 


