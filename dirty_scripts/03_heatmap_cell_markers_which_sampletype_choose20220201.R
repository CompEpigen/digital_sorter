library(dplyr)
library("viridis") 
library(ggplot2)
library(UpSetR)
# Set working directory ####
wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)
rawdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/patient_hlcma_addsplits/hlcma0001"
outdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/patient_hlcma_addsplits/hlcma0001"

# Set parameters ####
sample_type <- c("distal_normal","tan","luadc")

cell_in_interesed <- c("Club","Alveolar Epithelial Type 2")

#non_stratify_marker_selection_table####
## GO ####
c <- 2L
for(c in 1:length(cell_in_interesed)){
  cell_type <- cell_in_interesed[c]
  print(cell_type)
  
  
  markers_chose <- c()
  neg_markers_chose <- c()
  
  s <- 3L 
  for( s in 1:length(sample_type)){
    disease <- sample_type[s]
    print(disease)
   
    
    ## read files 
    table <- readRDS(paste0(outdir,"/non_stratify_marker_selection_table/",sample_type[s],"_marker_gene_table_subset_arranged.rds"))
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
  heatmap <- as.data.frame(heatmap)
  #row.names(heatmap) <- heatmap$markers
  #heatmap$markers <- NULL
  row.names.remove <- c("No markers selected")
  heatmap <-heatmap[!(heatmap$markers %in% row.names.remove),]
  club_colname <- c("markers Club",paste(sample_type,"Club") )
  AT2_colname <- c("markers AT2",paste(sample_type,"AT2"))
  # Club
  if(c==1){
    heatmap_club <- heatmap
    colnames(heatmap_club) <- club_colname
  }
  
  # AT2
  if(c==2){
    heatmap_AT2 <- heatmap
    colnames(heatmap_AT2) <- AT2_colname
  }
} 
  heatmap <- merge(heatmap_club,heatmap_AT2,
                   by.x = "markers Club",
                   by.y = "markers AT2",all =T)
  
  row.names(heatmap) <- heatmap$"markers Club"
  

 
  heatmap$"markers Club" <- NULL
  heatmap[is.na(heatmap)] <- 0
 
  df <- data.frame(x = colnames(heatmap))
  annotation <- as.data.frame(stringr::str_split_fixed(df$x, " ", 2))
  
  row.names(annotation) <- colnames(heatmap)
  colnames(annotation) <- c("Sample Type","Cell Type")
  Var1 = c("navy", "pink2", "lightblue")
  names(Var1) = sample_type
  Var3 <- c("moccasin", "firebrick")
  names(Var3) = c("Club","AT2")
  
  
  ann_colors = list( "Sample Type" = Var1,"Cell Type"= Var3)
  
  pdf(paste0(outdir,"/non_stratify_marker_selection_table/Heatmap_GOmarker_chose_clusterRow_AT2_Club_compare.pdf"), width=6, height=8)
  
  a <- pheatmap::pheatmap(heatmap,
                          main = "GO Markers, our own dataset",
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

## SurfaceGenie ##############
  #c <- 2L
  for(c in 1:length(cell_in_interesed)){
    cell_type <- cell_in_interesed[c]
    print(cell_type)
    
    
    markers_chose <- c()
    neg_markers_chose <- c()
    
    s <- 3L 
    for( s in 1:length(sample_type)){
      disease <- sample_type[s]
      print(disease)
      
      
      ## read files 
      table <- readRDS(paste0(outdir,"/non_stratify_marker_selection_table/",sample_type[s],"_SurfaceGenie_marker_gene_table_subset_arranged.rds"))
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
    heatmap <- as.data.frame(heatmap)
    #row.names(heatmap) <- heatmap$markers
    #heatmap$markers <- NULL
    row.names.remove <- c("No markers selected")
    heatmap <-heatmap[!(heatmap$markers %in% row.names.remove),]
    club_colname <- c("markers Club",paste(sample_type,"Club") )
    AT2_colname <- c("markers AT2",paste(sample_type,"AT2"))
    # Club
    if(c==1){
      heatmap_club <- heatmap
      colnames(heatmap_club) <- club_colname
    }
    
    # AT2
    if(c==2){
      heatmap_AT2 <- heatmap
      colnames(heatmap_AT2) <- AT2_colname
    }
  } 
  Sheatmap <- merge(heatmap_club,heatmap_AT2,
                   by.x = "markers Club",
                   by.y = "markers AT2",all =T) #78
  
  row.names(Sheatmap) <- Sheatmap$"markers Club"
  
  
  
  Sheatmap$"markers Club" <- NULL
  Sheatmap[is.na(Sheatmap)] <- 0
  
  df <- data.frame(x = colnames(Sheatmap))
  annotation <- as.data.frame(stringr::str_split_fixed(df$x, " ", 2))
  
  row.names(annotation) <- colnames(Sheatmap)
  colnames(annotation) <- c("Sample Type","Cell Type")
  Var1 = c("navy", "pink2", "lightblue")
  names(Var1) = sample_type
  Var3 <- c("moccasin", "firebrick")
  names(Var3) = c("Club","AT2")
  
  
  ann_colors = list( "Sample Type" = Var1,"Cell Type"= Var3)
  
  pdf(paste0(outdir,"/non_stratify_marker_selection_table/Heatmap_SurfaceGenie_marker_chose_clusterRow_AT2_Club_compare.pdf"), width=7, height=12)
  
  a <- pheatmap::pheatmap(Sheatmap,
                          main = "SurfaceGenie Markers, our own dataset",
                          
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
  

## SurfaceGenie Vs GO #######
intersect2 <- intersect(row.names(heatmap),row.names(Sheatmap))
new_heatmap <- Sheatmap[!(row.names(Sheatmap)%in%intersect2 ),]
pdf(paste0(outdir,"/non_stratify_marker_selection_table/Heatmap_NewS_marker_chose_clusterRow_AT2_Club_compare.pdf"), width=6, height=9)

a <- pheatmap::pheatmap(new_heatmap,
                        main = "New SurfaceGenie Markers, our own dataset",
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





# Level 2 ####
## GO ####
#c <- 1L
for(c in 1:length(cell_in_interesed)){
  cell_type <- cell_in_interesed[c]
  print(cell_type)
  
  
  markers_chose <- c()
  neg_markers_chose <- c()
  s <- 1L 
  for( s in 1:length(sample_type)){
    disease <- sample_type[s]
    print(disease)
    
    
    ## read files 
    table <- readRDS(paste0(outdir,"/level2_marker_selection_table/",sample_type[s],"_marker_gene_table_subset_arranged.rds"))
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
  heatmap <- as.data.frame(heatmap)
  #row.names(heatmap) <- heatmap$markers
  #heatmap$markers <- NULL
  row.names.remove <- c("No markers selected")
  heatmap <-heatmap[!(heatmap$markers %in% row.names.remove),]
  club_colname <- c("markers Club",paste(sample_type,"Club") )
  AT2_colname <- c("markers AT2",paste(sample_type,"AT2"))
  # Club
  if(c==1){
    heatmap_club <- heatmap
    colnames(heatmap_club) <- club_colname
  }
  
  # AT2
  if(c==2){
    heatmap_AT2 <- heatmap
    colnames(heatmap_AT2) <- AT2_colname
  }
} 

heatmap <- merge(heatmap_club,heatmap_AT2,
                 by.x = "markers Club",
                 by.y = "markers AT2",all =T)

row.names(heatmap) <- heatmap$"markers Club"



heatmap$"markers Club" <- NULL
heatmap[is.na(heatmap)] <- 0

df <- data.frame(x = colnames(heatmap))
annotation <- as.data.frame(stringr::str_split_fixed(df$x, " ", 2))

row.names(annotation) <- colnames(heatmap)
colnames(annotation) <- c("Sample Type","Cell Type")
Var1 = c("navy", "pink2", "lightblue")
names(Var1) = sample_type
Var3 <- c("moccasin", "firebrick")
names(Var3) = c("Club","AT2")


ann_colors = list( "Sample Type" = Var1,"Cell Type"= Var3)

pdf(paste0(outdir,"/level2_marker_selection_table/Heatmap_GOmarker_level2_clusterRow_AT2_Club_compare.pdf"), width=6, height=7)

a <- pheatmap::pheatmap(heatmap,
                        main = "GO Markers, our own dataset (level2)",
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

## SurfaceGenie ##############
c <- 2L
for(c in 1:length(cell_in_interesed)){
  cell_type <- cell_in_interesed[c]
  print(cell_type)
  
  
  markers_chose <- c()
  neg_markers_chose <- c()
  #s <- 1L 
  for( s in 1:length(sample_type)){
    disease <- sample_type[s]
    print(disease)
    
    
    #read files 
    table <- readRDS(paste0(outdir,"/level2_marker_selection_table/",sample_type[s],"_SurfaceGenie_marker_gene_table_subset_arranged.rds"))
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
  heatmap <- as.data.frame(heatmap)
  #row.names(heatmap) <- heatmap$markers
  #heatmap$markers <- NULL
  row.names.remove <- c("No markers selected")
  heatmap <-heatmap[!(heatmap$markers %in% row.names.remove),]
  club_colname <- c("markers Club",paste(sample_type,"Club") )
  AT2_colname <- c("markers AT2",paste(sample_type,"AT2"))
  # Club
  if(c==1){
    heatmap_club <- heatmap
    colnames(heatmap_club) <- club_colname
  }
  
  # AT2
  if(c==2){
    heatmap_AT2 <- heatmap
    colnames(heatmap_AT2) <- AT2_colname
  }
} 

Sheatmap <- merge(heatmap_club,heatmap_AT2,
                  by.x = "markers Club",
                  by.y = "markers AT2",all =T)

row.names(Sheatmap) <- Sheatmap$"markers Club"



Sheatmap$"markers Club" <- NULL
Sheatmap[is.na(Sheatmap)] <- 0

df <- data.frame(x = colnames(Sheatmap))
annotation <- as.data.frame(stringr::str_split_fixed(df$x, " ", 2))

row.names(annotation) <- colnames(Sheatmap)
colnames(annotation) <- c("Sample Type","Cell Type")
Var1 = c("navy", "pink2", "lightblue")
names(Var1) = sample_type
Var3 <- c("moccasin", "firebrick")
names(Var3) = c("Club","AT2")


ann_colors = list( "Sample Type" = Var1,"Cell Type"= Var3)

pdf(paste0(outdir,"/level2_marker_selection_table/Heatmap_SurfaceGenie_marker_level2_clusterRow_AT2_Club_compare.pdf"), width=7, height=10)

a <- pheatmap::pheatmap(Sheatmap,
                        main = "SurfaceGenie Markers, our own dataset (level2)",
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



## SurfaceGenie Vs GO #######
intersect2 <- intersect(row.names(heatmap),row.names(Sheatmap))
new_heatmap <- Sheatmap[!(row.names(Sheatmap)%in%intersect2 ),]
pdf(paste0(outdir,"/level2_marker_selection_table/Heatmap_NewS_marker_level2_clusterRow_AT2_Club_compare.pdf"), width=7, height=6)

a <- pheatmap::pheatmap(new_heatmap,
                        main = "New SurfaceGenie Markers, our own dataset (level 2)",
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




# expression heatmap####
features <- rownames(Sheatmap) #using non-stratified

heatmap_data <- data.frame("gene"=features)
empty <- data.frame("gene"=features)

ob <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/patient_hlcma/patient_hlcma_addsplits.rds")

split <- as.character(unique(sort(ob@meta.data$disease)))
ob_split <- SplitObject(ob, split.by = "disease")
length(split)

#k <- 1L
for(k in 1:length(split)){
  print(split[k])      
  ob_split2 = ob_split[[split[k]]]
  ob_split2 <- SplitObject(ob_split2, split.by = "split1")
  
  ob2 = ob_split2$`PTPRC-`
  ob22 <- SplitObject(ob2, split.by = "split2")
  ob_split_PTPRC0 = ob22$`EPCAM+`
  
  a <- DotPlot(ob_split_PTPRC0, features = features, group.by = "annotation.l2") + 
    RotatedAxis()+  
    theme(plot.title = element_text(hjust = 0.5)) 
  data <- a[["data"]]
  data <- data[,c(1,3,4)]
  cell_type <- unique(data$id)
  data2 <- empty
  #c <- 2L
  for(c in seq(length(cell_type))){
    print(cell_type[c])
    data_cell <- data %>% filter(id == cell_type[c])
    data_cell$id <- NULL
    colnames(data_cell) <- c(paste0(cell_type[c],":",split[k]),"gene")
    data2 <- merge(data2,data_cell,by="gene",all=T) 
  }
  
  
  
  heatmap_data <- merge(heatmap_data,data2,by="gene",all=T)
  rm(ob_split_PTPRC0)
  rm(ob2) 
  rm(ob22)
}    
  
heatmap_data[is.na(heatmap_data)] <- 0
heatmap_data <- as.data.frame(heatmap_data)
row.names(heatmap_data) <- heatmap_data$"gene"
heatmap_data$"gene"<- NULL
saveRDS(heatmap_data,paste0(outdir,"/20220201_Heatmap_expression_data_all_level2_SurfaceGenie_genes.rds"))

## can start here!!!!!!!#############
heatmap_data <- readRDS(paste0(outdir,"/20220201_Heatmap_expression_data_all_level2_SurfaceGenie_genes.rds"))
df <- data.frame(x = colnames(heatmap_data))
annotation <- as.data.frame(stringr::str_split_fixed(df$x, ":", 2))
row.names(annotation) <- colnames(heatmap_data)
colnames(annotation) <- c("Cell Type","Sample Type")
annotation$`Cell Type` <- gsub("Signaling Alveolar Epithelial Type 2","Alveolar Epithelial Type 2",annotation$`Cell Type`)

Var2 = c("navy", "lightblue","pink2")
names(Var2) = sample_type
Var3 <- c("hotpink", "firebrick","green2","gold","moccasin","darkgreen","black",
          "steelblue1","gray","navy","limegreen","orange","purple1")
names(Var3) = unique(sort(annotation$"Cell Type"))

ann_colors = list("Sample Type" = Var2,"Cell Type"= Var3)


### non- pearson ##########
cell_order <- c("Signaling Alveolar Epithelial Type 2","Alveolar Epithelial Type 2" ,
                "Club","Alveolar Epithelial Type 1",
                "Mucous" ,"Goblet","Basal","Differentiating Basal","Proliferating Basal",
                "Proximal Basal", 
                "Ciliated" ,"Proximal Ciliated","Ionocyte","Neuroendocrine"  ,"Serous")
heatmap_sub <- heatmap_data %>% dplyr:: select(starts_with(cell_order))

pdf(paste0(outdir,"/20220201_Heatmap_scale_experssion_allLevel2_clusterRow_euclidean.pdf"), width=12, height=14)

a <- pheatmap::pheatmap(heatmap_sub,
                        color = colorRampPalette(c("blue4","white", "red2"))(50),
                        main = "SurfaceGenie Markers chosen for Club and AT2 cells",
                        annotation = annotation,
                        border_color = "black", 
                        show_colnames=F,
                        cellheight=10,
                        scale = "row",
                        annotation_colors = ann_colors, 
                        clustering_distance_rows = "euclidean", 
                        cluster_rows = T,
                        cluster_cols = F)
print(a)
dev.off()

pdf(paste0(outdir,"/20220201_Heatmap_scale_experssion_allLevel2_clusterRow_Col_euclidean.pdf"), width=12, height=14)

a <- pheatmap::pheatmap(heatmap_sub,
                        color = colorRampPalette(c("blue4","white", "red2"))(50),
                        main = "SurfaceGenie Markers chosen for Club and AT2 cells",
                        annotation = annotation,
                        border_color = "black", 
                        show_colnames=F,
                        cellheight=10,
                        scale = "row",
                        annotation_colors = ann_colors, 
                        clustering_distance_rows = "euclidean", 
                        cluster_rows = T,
                        cluster_cols = T)
print(a)
dev.off()

