library(dplyr)
library("viridis") 
library(ggplot2)
library(UpSetR)
library(Seurat)
library(SeuratObject)
# Set working directory ####
wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)
rawdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/patient_hlcma_addsplits/hlcma0001"
outdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/patient_hlcma_addsplits/hlcma0001"
# load colors
library(randomcoloR)
#n <- 50
#palette <- distinctColorPalette(n)
#saveRDS(palette,"/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/50color_palette.RDS")
palette <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/50color_palette.RDS")
# Set parameters ####
ob <- readRDS(paste0("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/patient_hlcma/patient_hlcma_tumor_added_addsplits_fixedAnno.rds"))
sampledata <- ob@meta.data

sample_type <- unique(sort(sampledata$disease))
sample_type <- sample_type[c(2,1,4,3)]
#"general_normal" "distal_normal"  "tan"            "luadc"        


ob_split <- SplitObject(ob, split.by = "split1")
ob_split1<-ob_split[["PTPRC+"]]
ob_split2<-ob_split[["PTPRC-"]]
ob_split2 <- SplitObject(ob_split2, split.by = "split2")
ob_split2 <- ob_split2[["EPCAM+"]]
#non_stratify_marker_selection_table####
immune_cells <- unique(sort(ob_split1@meta.data[["annotation.l2"]]))
epithelial_cells <-unique(sort(ob_split2@meta.data[["annotation.l2"]]))

## SurfaceGenie ##############
### immune cells ####
  #c <- 7L
for(c in 1:length(immune_cells)){
  cell_type <- immune_cells[c]
  print(cell_type)
  
  
  markers_chose <- c()
  neg_markers_chose <- c()
  
  s <- 3L 
  for( s in 1:length(sample_type)){
    disease <- sample_type[s]
    print(disease)
    
    
    ## read files 
    table <- readRDS(paste0(outdir,"/non_stratify_marker_selection_table/add_background/",sample_type[s],"_table_add_background.rds"))
    table <- table %>% filter(table$annotation.l2 == cell_type)
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
  
  heatmap_colname <- c("markers",paste0(sample_type,":",immune_cells[c]) )
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


Var1 = c("navy", "royalblue", "lightblue", "pink2")
names(Var1) = sample_type

Var2 = palette[seq(length(immune_cells))]
names(Var2) = immune_cells
ann_colors = list( "Sample Type" = Var1,"Cell Type"= Var2)
  
  
  
  
  
df <- data.frame(x = colnames(Sheatmap))
annotation <- as.data.frame(stringr::str_split_fixed(df$x, ":", 2))
  
row.names(annotation) <- colnames(Sheatmap)
colnames(annotation) <- c("Sample Type","Cell Type")


  
pdf(paste0(outdir,"/non_stratify_marker_selection_table/Heatmap_Immune_cells_compare.pdf"), width=12, height=42)
  
a <- pheatmap::pheatmap(Sheatmap,
                        main = "SurfaceGenie Markers, our own dataset with tumor annotation",
                        
                        color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                        annotation = annotation,
                        legend_breaks= c(1,0,-1),
                        legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                        annotation_colors = ann_colors, 
                        cluster_rows = T,
                        cluster_cols = F,
                        show_colnames= F)
                         
  
  print(a)
  
  b <- pheatmap::pheatmap(Sheatmap,
                          main = "SurfaceGenie Markers, our own dataset with tumor annotation",
                          
                          color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                          annotation = annotation,
                          legend_breaks= c(1,0,-1),
                          legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                          annotation_colors = ann_colors, 
                          cluster_rows = T,
                          cluster_cols = T,
                          show_colnames= F)
  
  
  print(b)
  
  dev.off()
  


### epithelial cells ####
#c <- 7L
  for(c in 1:length(epithelial_cells)){
    cell_type <- epithelial_cells[c]
    print(cell_type)
    
    
    markers_chose <- c()
    neg_markers_chose <- c()
    
    s <- 3L 
    for( s in 1:length(sample_type)){
      disease <- sample_type[s]
      print(disease)
      
      
      ## read files 
      table <- readRDS(paste0(outdir,"/non_stratify_marker_selection_table/add_background/",sample_type[s],"_table_add_background.rds"))
      table <- table %>% filter(table$annotation.l2 == cell_type)
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
    
    heatmap_colname <- c("markers",paste0(sample_type,":",epithelial_cells[c]) )
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
  annotation$`Cell Type` <- gsub("Signaling Alveolar Epithelial Type 2","Alveolar Epithelial Type 2",annotation$`Cell Type`)
  
  Var1 = c("navy", "royalblue", "lightblue", "pink2")
  names(Var1) = sample_type
  
  Var2 =  c("hotpink", "firebrick","green2",
            "blue2","royalblue4","steelblue4","deepskyblue",
            "grey","gold","moccasin","darkgreen","black",
            "lightblue","orangered","navy","limegreen","turquoise","orange","purple1")
  names(Var2) =unique(sort(annotation$"Cell Type"))
  ann_colors = list( "Sample Type" = Var1,"Cell Type"= Var2)  
  
  cell_order <- c("C9_Unknown", "C1_Tumor","C4_Tumor","C5_Tumor","C15_Tumor",
                  "Signaling Alveolar Epithelial Type 2","Alveolar Epithelial Type 2" ,
                  "Club","Alveolar Epithelial Type 1",
                  "Mucous" ,"Goblet","Basal","Differentiating Basal","Proliferating Basal",
                  "Proximal Basal", 
                  "Ciliated" ,"Proximal Ciliated","Ionocyte","Neuroendocrine"  ,"Serous")
  Sheatmap <- Sheatmap %>% dplyr:: select(contains(cell_order))
  
pdf(paste0(outdir,"/non_stratify_marker_selection_table/Heatmap_epithelial_cells_compare.pdf"), width=16, height=40)
  
  a <- pheatmap::pheatmap(Sheatmap,
                          main = "SurfaceGenie Markers, our own dataset with tumor annotation",
                          
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
                          main = "SurfaceGenie Markers, our own dataset with tumor annotation",
                          
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
  
##special heatmap ####  
  # show only gene_interest
  #gene_interest <- c("SLC22A4","TACSTD2","AQP5","TSPAN8","ANXA1","LGALS3","MUC16","CD24","TSPAN1","BSG","ITGB4","PLAUR","CX3CL1","KRTCAP2","BST2","DCBLD2","AZGP1","ABCC3","PLAU","ADGRF1","MFGE8","CEACAM6","CEACAM5","SERINC2","EPCAM","MSLN","TNFRSF12A","GPRC5A","TM4SF1","MIF","KRT18","ATP1B1","CD63")
  
  tumor_marker <- c("STEAP4", "PLAUR", "DCBLD2","ITGB4",  "ITGA3",  "KRTCAP2","MUC16","CD24","HMGB1", "TSPAN1",
                    "ATP5IF1", "TMEM123", "CD59",  "SLC44A4","F11R","CEACAM5","CEACAM6","CD46","CD63","MSLN",
                    "SERINC2","TSPAN8","CLDN3","TMBIM6","AQP5","CRLF1","PDIA3","HSPD1","ENO1","P4HB",
                    "KRT18","SPINT2","ITGA6","PIEZO1","ABCC3","CX3CL1","AZGP1","CD274","BSG","BST2","MFGE8",
                    "HCAR2","ATRAID","CTSZ","TMED10","KRT10","TMBIM4","HMMR",
                    "LRP5","HACD3","TM4SF1"
                    )
  
 

  ##attention!!!!which Sheatmap using!
  Sheatmap2 <- Sheatmap%>% filter(rownames(Sheatmap)%in% tumor_marker) #######
  Sheatmap2 <- Sheatmap2[tumor_marker,]
  Sheatmap2 <- na.omit(Sheatmap2)
  pdf(paste0(outdir,"/non_stratify_marker_selection_table/Heatmap_TumorMarker_epithelial_cells.pdf"), width=12, height=8)
  
  a <- pheatmap::pheatmap(Sheatmap2,
                          main = "Potential Tumor Marker, our own dataset with tumor annotation",
                          
                          color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                          annotation = annotation,
                          legend_breaks= c(1,0,-1),
                          legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                          annotation_colors = ann_colors, 
                          cluster_rows = T,
                          cluster_cols = F,
                          show_colnames= F)
  
  
  print(a)
  
  dev.off()
  
  
  ### show only genes appears less than 10 times ####
  
  Sheatmap$rowsum <- rowSums(Sheatmap)
  Sheatmap3 <- Sheatmap%>% filter(abs(rowsum) < 10)
  Sheatmap3$rowsum <- NULL
  pdf(paste0(outdir,"/non_stratify_marker_selection_table/Heatmap_unique_epithelial_cells.pdf"), width=14, height=35)
  
  a <- pheatmap::pheatmap(Sheatmap3,
                          main = "SurfaceGenie Markers (not shared by so many cell types), our own dataset with tumor annotation",
                          
                          color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                          annotation = annotation,
                          legend_breaks= c(1,0,-1),
                          legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                          annotation_colors = ann_colors, 
                          cluster_rows = T,
                          cluster_cols = F,
                          show_colnames= T)
  
  
  print(a)
  
 dev.off()

 
 # Level 2 ####
## SurfaceGenie ##############
 ### epithelial cells ####
 #c <- 7L
 for(c in 1:length(epithelial_cells)){
   cell_type <- epithelial_cells[c]
   print(cell_type)
   
   
   markers_chose <- c()
   neg_markers_chose <- c()
   
   s <- 3L 
   for( s in 1:length(sample_type)){
     disease <- sample_type[s]
     print(disease)
     
     
     ## read files 
     table <- readRDS(paste0(outdir,"/level2_marker_selection_table/add_background/",sample_type[s],"_table_add_background.rds"))
     table <- table %>% filter(table$annotation.l2 == cell_type)
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
   
   heatmap_colname <- c("markers",paste0(sample_type,":",epithelial_cells[c]) )
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
 annotation$`Cell Type` <- gsub("Signaling Alveolar Epithelial Type 2","Alveolar Epithelial Type 2",annotation$`Cell Type`)
 
 Var1 = c("navy", "royalblue", "lightblue", "pink2")
 names(Var1) = sample_type
 
 Var2 =  c("hotpink", "firebrick","green2",
           "blue2","royalblue4","steelblue4","deepskyblue",
           "grey","gold","moccasin","darkgreen","black",
           "lightblue","orangered","navy","limegreen","turquoise","orange","purple1")
 names(Var2) =unique(sort(annotation$"Cell Type"))
 ann_colors = list( "Sample Type" = Var1,"Cell Type"= Var2)
 
 cell_order <- c("C9_Unknown", "C1_Tumor","C4_Tumor","C5_Tumor","C15_Tumor",
                 "Signaling Alveolar Epithelial Type 2","Alveolar Epithelial Type 2" ,
                 "Club","Alveolar Epithelial Type 1",
                 "Mucous" ,"Goblet","Basal","Differentiating Basal","Proliferating Basal",
                 "Proximal Basal", 
                 "Ciliated" ,"Proximal Ciliated","Ionocyte","Neuroendocrine"  ,"Serous")
 Sheatmap <- Sheatmap %>% dplyr:: select(contains(cell_order))
 
 
pdf(paste0(outdir,"/level2_marker_selection_table/Heatmap_level2_epithelial_cells.pdf"), width=16, height=38)
 
 a <- pheatmap::pheatmap(Sheatmap,
                         main = "SurfaceGenie Markers, our own dataset with tumor annotation (Level2)",
                         
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
                         main = "SurfaceGenie Markers, our own dataset with tumor annotation (Level2)",
                         
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
 
 # show only gene_interest
 Sheatmap2 <- Sheatmap%>% filter(rownames(Sheatmap)%in% tumor_marker)
 Sheatmap2 <- Sheatmap2[tumor_marker,]
 Sheatmap2 <- na.omit(Sheatmap2)

pdf(paste0(outdir,"/level2_marker_selection_table/Heatmap_TumorMarker_epithelial_cells.pdf"), width=14, height=12)
 
 a <- pheatmap::pheatmap(Sheatmap2,
                         main = "SurfaceGenie Markers, our own dataset with tumor annotation (Level2)",
                         
                         color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                         annotation = annotation,
                         legend_breaks= c(1,0,-1),
                         legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                         annotation_colors = ann_colors, 
                         cluster_rows = T,
                         cluster_cols = F,
                         show_colnames= T)
 
 
 print(a)
 
 dev.off()
 
 # show only genes appears less than 10 times
 
 Sheatmap$rowsum <- rowSums(Sheatmap)
 Sheatmap3 <- Sheatmap%>% filter(abs(rowsum) < 10)
 Sheatmap3$rowsum <- NULL
 
pdf(paste0(outdir,"/level2_marker_selection_table/Heatmap_unique_epithelial_cells.pdf"), width=14, height=35)
 
 a <- pheatmap::pheatmap(Sheatmap3,
                         main = "SurfaceGenie Markers, our own dataset with tumor annotation (Level2)",
                         
                         color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                         annotation = annotation,
                         legend_breaks= c(1,0,-1),
                         legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                         annotation_colors = ann_colors, 
                         cluster_rows = T,
                         cluster_cols = F,
                         show_colnames= T)
 
 
 print(a)
 
 dev.off()

 # CD45- ####
 ## SurfaceGenie ##############
 ### epithelial cells ####
 #c <- 7L
 for(c in 1:length(epithelial_cells)){
   cell_type <- epithelial_cells[c]
   print(cell_type)
   
   
   markers_chose <- c()
   neg_markers_chose <- c()
   
   s <- 3L 
   for( s in 1:length(sample_type)){
     disease <- sample_type[s]
     print(disease)
     
     
     ## read files 
     table <- readRDS(paste0(outdir,"/CD45-_marker_selection_table/add_background/",sample_type[s],"_table_add_background.rds"))
     table <- table %>% filter(table$annotation.l2 == cell_type)
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
   
   heatmap_colname <- c("markers",paste0(sample_type,":",epithelial_cells[c]) )
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
 annotation$`Cell Type` <- gsub("Signaling Alveolar Epithelial Type 2","Alveolar Epithelial Type 2",annotation$`Cell Type`)
 
 Var1 = c("navy", "royalblue", "lightblue", "pink2")
 names(Var1) = sample_type
 
 Var2 =  c("hotpink", "firebrick","green2",
           "blue2","royalblue4","steelblue4","deepskyblue",
           "grey","gold","moccasin","darkgreen","black",
           "lightblue","orangered","navy","limegreen","turquoise","orange","purple1")
 names(Var2) =unique(sort(annotation$"Cell Type"))
 ann_colors = list( "Sample Type" = Var1,"Cell Type"= Var2)
 
 cell_order <- c("C9_Unknown", "C1_Tumor","C4_Tumor","C5_Tumor","C15_Tumor",
                 "Signaling Alveolar Epithelial Type 2","Alveolar Epithelial Type 2" ,
                 "Club","Alveolar Epithelial Type 1",
                 "Mucous" ,"Goblet","Basal","Differentiating Basal","Proliferating Basal",
                 "Proximal Basal", 
                 "Ciliated" ,"Proximal Ciliated","Ionocyte","Neuroendocrine"  ,"Serous")
 Sheatmap <- Sheatmap %>% dplyr:: select(contains(cell_order))
 
 pdf(paste0(outdir,"/CD45-_marker_selection_table/Heatmap_CD45-_epithelial_cells.pdf"), width=16, height=40)
 
 a <- pheatmap::pheatmap(Sheatmap,
                         main = "SurfaceGenie Markers, our own dataset with tumor annotation (CD45-)",
                         
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
                         main = "SurfaceGenie Markers, our own dataset with tumor annotation (CD45-)",
                         
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

 # show only gene_interest
 Sheatmap2 <- Sheatmap%>% filter(rownames(Sheatmap)%in% tumor_marker)
 Sheatmap2 <- Sheatmap2[tumor_marker,]
 Sheatmap2 <- na.omit(Sheatmap2)
pdf(paste0(outdir,"/CD45-_marker_selection_table/Heatmap_TumorMarker_epithelial_cells.pdf"), width=14, height=10)
 
 a <- pheatmap::pheatmap(Sheatmap2,
                         main = "SurfaceGenie Markers, our own dataset with tumor annotation (CD45-)",
                         
                         color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                         annotation = annotation,
                         legend_breaks= c(1,0,-1),
                         legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                         annotation_colors = ann_colors, 
                         cluster_rows = F,
                         cluster_cols = F,
                         show_colnames= T)
 
 
 print(a)
 
 dev.off()
 # show only genes appears less than 10 times
 
 Sheatmap$rowsum <- rowSums(Sheatmap)
 Sheatmap3 <- Sheatmap%>% filter(abs(rowsum) < 5)
 Sheatmap3$rowsum <- NULL
 
 pdf(paste0(outdir,"/CD45-_marker_selection_table/Heatmap_unique5_epithelial_cells.pdf"), width=14, height=30)
 
 a <- pheatmap::pheatmap(Sheatmap3,
                         main = "SurfaceGenie Markers, our own dataset with tumor annotation (CD45-)",
                         
                         color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                         annotation = annotation,
                         legend_breaks= c(1,0,-1),
                         legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                         annotation_colors = ann_colors, 
                         cluster_rows = T,
                         cluster_cols = F,
                         show_colnames= T)
 
 
 print(a)
 
 dev.off()
 
