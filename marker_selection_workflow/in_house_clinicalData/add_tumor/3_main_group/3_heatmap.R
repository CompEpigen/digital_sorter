library(dplyr)
library("viridis") 
library(ggplot2)
library(UpSetR)
library(Seurat)
library(SeuratObject)
# Set working directory ####
wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)
rawdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/patient_hlcma"

outdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/patient_hlcma_addsplits/hlcma0001/3groups_level2"
# load colors
library(randomcoloR)
#n <- 50
#palette <- distinctColorPalette(n)
#saveRDS(palette,"/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/50color_palette.RDS")
palette <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/50color_palette.RDS")
# Set parameters ####

group1 <- c("Tumor","Alveolar_Epithelial")
group2 <- c("Tumor","Airway_Epithelial")

sample_type <- c("distal_normal","tan","luadc")
#non_stratify_marker_selection_table####


## SurfaceGenie ##############

### epithelial cells ####
#c <- 7L
  for(c in 1:length(group1)){
    cell_type <- group1[c]
    print(cell_type)
    
    
    markers_chose <- c()
    neg_markers_chose <- c()
    
    
    for( s in 1:length(sample_type)){
      disease <- sample_type[s]
      print(disease)
      
      
      ## read files 
      table <- readRDS(paste0(outdir,"/add_background/",sample_type[s],"_marker_gene_table_Tumor_Vs_Alveolar_background.rds"))
      table <- table %>% filter(table$annotation.l3 == cell_type)
      pos_table <- table %>% filter(table$avg_log2FC >0)
      neg_table <- table %>% filter(table$avg_log2FC <0)
      markers <- pos_table$gene
      markers <- unique(sort(markers))
      
      neg_markers <- neg_table$gene
      neg_markers <- unique(sort(neg_markers))
      if(length(markers)==0){
        markers <- "No markers selected"
      }
      if(length(neg_markers)==0){
        neg_markers <- "No markers selected"
      }
      markers2 <- as.data.frame(cbind (markers=markers,group = rep(sample_type[s],length(markers))))
      markers3 <- as.data.frame(cbind (markers=neg_markers,group = rep(sample_type[s],length(neg_markers))))
      markers_chose <- rbind(markers_chose,markers2)
      neg_markers_chose <- rbind(neg_markers_chose,markers3)
      
    }
    
    markers_chose$value <- rep(1,dim(markers_chose)[1])
    neg_markers_chose$value <- rep(-1,dim(neg_markers_chose)[1])
    heatmap_pos <- as.data.frame(tidyr::pivot_wider(markers_chose, names_from = group, values_from = value))
    heatmap_neg <- as.data.frame(tidyr::pivot_wider(neg_markers_chose, names_from = group, values_from = value))
    
    heatmap <- rbind(heatmap_pos,heatmap_neg)
    heatmap[is.na(heatmap)] <- 0
    heatmap <- aggregate(. ~ markers, heatmap,sum)
    
    heatmap <- as.data.frame(heatmap)
    #row.names(heatmap) <- heatmap$markers
    #heatmap$markers <- NULL
    row.names.remove <- c("No markers selected")
    heatmap <-heatmap[!(heatmap$markers %in% row.names.remove),]
    
    heatmap_colname <- c("markers",paste0(sample_type,":",group1[c]) )
    # Club
    if(c==1){
      Sheatmap <- heatmap
      colnames(Sheatmap) <- heatmap_colname
    }else if(c!=1){
      colnames(heatmap) <- heatmap_colname
      Sheatmap <- merge(Sheatmap,heatmap,
                        by = "markers", all =T) 
      
    }
  }  
  row.names(Sheatmap) <- Sheatmap$"markers"
  Sheatmap$"markers" <- NULL
  Sheatmap[is.na(Sheatmap)] <- 0

 
  df <- data.frame(x = colnames(Sheatmap))
  annotation <- as.data.frame(stringr::str_split_fixed(df$x, ":", 2))
  
  row.names(annotation) <- colnames(Sheatmap)
  colnames(annotation) <- c("Sample Type","Cell Type")
 
  Var1 = c("royalblue", "lightblue", "pink2")
  names(Var1) = sample_type
  
  Var2 =  c("orange","darkblue")
  names(Var2) =unique(sort(annotation$"Cell Type"))
  ann_colors = list( "Sample Type" = Var1,"Cell Type"= Var2)  
  
  
  
  
pdf(paste0(outdir,"/Heatmap_Tumor_Alveolar_compare.pdf"), width=10, height=10)
  
  a <- pheatmap::pheatmap(Sheatmap,
                          main = "Pan-Cancer's Cluster Markers, our own dataset with tumor annotation",
                          
                          color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                          annotation = annotation,
                          legend_breaks= c(1,0,-1),
                          legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                          annotation_colors = ann_colors, 
                          cluster_rows = T,
                          cluster_cols = F,
                          show_colnames= T)
  
  
  print(a)
  
  b <- pheatmap::pheatmap(Sheatmap,
                          main = "Pan-Cancer's Cluster Markers, our own dataset with tumor annotation",
                          
                          color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                          annotation = annotation,
                          legend_breaks= c(1,0,-1),
                          legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                          annotation_colors = ann_colors, 
                          cluster_rows = T,
                          cluster_cols = T,
                          show_colnames= T)
  
  
  print(b)
  
  dev.off() 
  
  # group2 ####
  for(c in 1:length(group2)){
    cell_type <- group2[c]
    print(cell_type)
    
    
    markers_chose <- c()
    neg_markers_chose <- c()
    
    
    for( s in 1:length(sample_type)){
      disease <- sample_type[s]
      print(disease)
      
      
      ## read files 
      table <- readRDS(paste0(outdir,"/add_background/",sample_type[s],"_marker_gene_table_Tumor_Vs_Airway_background.rds"))
      table <- table %>% filter(table$annotation.l3 == cell_type)
      pos_table <- table %>% filter(table$avg_log2FC >0)
      neg_table <- table %>% filter(table$avg_log2FC <0)
      markers <- pos_table$gene
      markers <- unique(sort(markers))
      
      neg_markers <- neg_table$gene
      neg_markers <- unique(sort(neg_markers))
      if(length(markers)==0){
        markers <- "No markers selected"
      }
      if(length(neg_markers)==0){
        neg_markers <- "No markers selected"
      }
      markers2 <- as.data.frame(cbind (markers=markers,group = rep(sample_type[s],length(markers))))
      markers3 <- as.data.frame(cbind (markers=neg_markers,group = rep(sample_type[s],length(neg_markers))))
      markers_chose <- rbind(markers_chose,markers2)
      neg_markers_chose <- rbind(neg_markers_chose,markers3)
      
    }
    
    markers_chose$value <- rep(1,dim(markers_chose)[1])
    neg_markers_chose$value <- rep(-1,dim(neg_markers_chose)[1])
    heatmap_pos <- as.data.frame(tidyr::pivot_wider(markers_chose, names_from = group, values_from = value))
    heatmap_neg <- as.data.frame(tidyr::pivot_wider(neg_markers_chose, names_from = group, values_from = value))
    
    heatmap <- rbind(heatmap_pos,heatmap_neg)
    heatmap[is.na(heatmap)] <- 0
    heatmap <- aggregate(. ~ markers, heatmap,sum)
    
    heatmap <- as.data.frame(heatmap)
    #row.names(heatmap) <- heatmap$markers
    #heatmap$markers <- NULL
    row.names.remove <- c("No markers selected")
    heatmap <-heatmap[!(heatmap$markers %in% row.names.remove),]
    
    heatmap_colname <- c("markers",paste0(sample_type,":",group2[c]) )
    # Club
    if(c==1){
      Sheatmap <- heatmap
      colnames(Sheatmap) <- heatmap_colname
    }else if(c!=1){
      colnames(heatmap) <- heatmap_colname
      Sheatmap <- merge(Sheatmap,heatmap,
                        by = "markers", all =T) 
      
    }
  }  
  row.names(Sheatmap) <- Sheatmap$"markers"
  Sheatmap$"markers" <- NULL
  Sheatmap[is.na(Sheatmap)] <- 0
  
  
  df <- data.frame(x = colnames(Sheatmap))
  annotation <- as.data.frame(stringr::str_split_fixed(df$x, ":", 2))
  
  row.names(annotation) <- colnames(Sheatmap)
  colnames(annotation) <- c("Sample Type","Cell Type")
  
  Var1 = c("royalblue", "lightblue", "pink2")
  names(Var1) = sample_type
  
  Var2 =  c("greenyellow","darkblue")
  names(Var2) =unique(sort(annotation$"Cell Type"))
  ann_colors = list( "Sample Type" = Var1,"Cell Type"= Var2)  
  
  
  
  
  pdf(paste0(outdir,"/Heatmap_Tumor_Airway_compare.pdf"), width=10, height=8)
  
  a <- pheatmap::pheatmap(Sheatmap,
                          main = "Pan-Cancer's Cluster Markers, our own dataset with tumor annotation",
                          
                          color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                          annotation = annotation,
                          legend_breaks= c(1,0,-1),
                          legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                          annotation_colors = ann_colors, 
                          cluster_rows = T,
                          cluster_cols = F,
                          show_colnames= T)
  
  
  print(a)
  
  b <- pheatmap::pheatmap(Sheatmap,
                          main = "Pan-Cancer's Cluster Markers, our own dataset with tumor annotation",
                          
                          color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                          annotation = annotation,
                          legend_breaks= c(1,0,-1),
                          legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                          annotation_colors = ann_colors, 
                          cluster_rows = T,
                          cluster_cols = T,
                          show_colnames= T)
  
  
  print(b)
  
  dev.off() 
  