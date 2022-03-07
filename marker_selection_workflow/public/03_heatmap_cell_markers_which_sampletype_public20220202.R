library(dplyr)
library("viridis") 
library(ggplot2)
library(UpSetR)
# Set working directory ####
wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)
rawdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/lung_mapped_cellxgene_fixed_addsplits"
outdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/lung_mapped_cellxgene_fixed_addsplits"

# Load functions
GetfileNames <- function(fileDir, pattern = ".rds"){
  filenames <- list.files(fileDir, pattern= pattern, full.names=TRUE)
  filename <- sub(paste0(".*",fileDir,"/"), "", filenames)
  filename <- sub(pattern, "", filename)  
  return(filename)
}

# Set parameters ####
list_sample_cohort <- readRDS(paste0(rawdir,"/list_sample_cohort.rds"))
sample_type <- c("general_normal",
                 "adjacent_normal","blood","luad","lusc",
                 "nsclc_primary","nsclc_meta","LymphNode_normal","LymphNode_meta","pleural_fluids")

cell_in_interesed <- c("Club","Alveolar Epithelial Type 2")

#non_stratify_marker_selection_table####
## GO ####
c <- 1L
for(c in 1:length(cell_in_interesed)){
  cell_type <- cell_in_interesed[c]
  print(cell_type)
  
  
  markers_chose <- c()
  neg_markers_chose <- c()
  
  #s <- 5L 
  for( s in 1:length(sample_type)){
    disease <- sample_type[s]
    print(disease)
    datasets <- list_sample_cohort[[disease]]
    
    #d <- 1L
    for(d in 1:length(datasets)){
      print(datasets[d])
      ## read files ####
      table <- readRDS(paste0(outdir,"/",datasets[d],"/",sample_type[s],"_marker_gene_table_subset_arranged.rds")) #_SurfaceGenie
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
      group <- paste(datasets[d],sample_type[s],sep = ":")
      markers2 <- as.data.frame(cbind (markers=markers,group = rep(group,length(markers))))
      markers3 <- as.data.frame(cbind (markers=neg_markers,group = rep(group,length(neg_markers))))
      markers_chose <- rbind(markers_chose,markers2)
      neg_markers_chose <- rbind(neg_markers_chose,markers3)
    }
   
    
  }

  markers_chose$value <- rep(1,dim(markers_chose)[1])
  neg_markers_chose$value <- rep(-1,dim(neg_markers_chose)[1])
  heatmap_pos <- as.data.frame(tidyr::pivot_wider(markers_chose, names_from = group, values_from = value))
  heatmap_neg <- as.data.frame(tidyr::pivot_wider(neg_markers_chose, names_from = group, values_from = value))
  
  heatmap <- rbind(heatmap_pos,heatmap_neg)
  
  heatmap[is.na(heatmap)] <- 0
  heatmap <- as.data.frame(heatmap)
  col_order <- c("markers","travaglini_2020:general_normal","bischoff_2021:adjacent_normal","kim_2020:adjacent_normal","lambrechts_2018:adjacent_normal", "song_2019:adjacent_normal",
                 "bischoff_2021:luad","lambrechts_2018:luad","wu_2021:luad",
                 "lambrechts_2018:lusc" ,"wu_2021:lusc",
                 "kim_2020:nsclc_primary","lambrechts_2018:nsclc_primary","wu_2021:nsclc_primary","zilionis_2019:nsclc_primary" ,"song_2019:nsclc_primary",
                 "kim_2020:nsclc_meta",
                 "kim_2020:pleural_fluids","kim_2020:LymphNode_normal", "kim_2020:LymphNode_meta",
                 "zilionis_2019:blood")
  heatmap <-aggregate(heatmap[,-1], by=list(markers=heatmap$markers), FUN=sum)
  heatmap <-heatmap[,col_order]
  row.names.remove <- c("No markers selected")
  heatmap <-heatmap[!(heatmap$markers %in% row.names.remove),]
  club_colname <- paste(col_order,"Club")
  AT2_colname <- paste(col_order,"AT2")
  
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
  annotation <- as.data.frame(stringr::str_split_fixed(df$x, ":", 2))
  annotation2 <- as.data.frame(stringr::str_split_fixed(annotation$V2, " ", 2))
  annotation <- cbind(annotation,annotation2)
  annotation <- annotation[,c(1,3,4)]
  row.names(annotation) <- colnames(heatmap)
  colnames(annotation) <- c("Cohort","Sample Type","Cell Type")
  Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
           "#661100", "#0072B2", "#D55E00")
  names(Var1) = unique(sort(annotation$Cohort))
  Var2 = c("lightblue", "darkred","navy","pink2","purple",
           "gold","orange","greenyellow","darkgreen","royalblue")
  names(Var2) = unique(sort(annotation$"Sample Type"))
  Var3 <- c("moccasin", "firebrick")
  names(Var3) = c("Club","AT2")
  
  ann_colors = list(Cohort = Var1, "Sample Type" = Var2,"Cell Type"= Var3)
  
  pdf(paste0(outdir,"/20220202_Heatmap_GOmarker_chose_clusterRow_AT2_Club_compare.pdf"), width=12, height=20)
  
  a <- pheatmap::pheatmap(heatmap,
                          main = "GO Markers, public datasets",
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
    
    s <- 3L 
    for( s in 1:length(sample_type)){
      disease <- sample_type[s]
      print(disease)
      datasets <- list_sample_cohort[[disease]]
      
      #d <- 1L
      for(d in 1:length(datasets)){
        print(datasets[d])
        ## read files ####
        table <- readRDS(paste0(outdir,"/",datasets[d],"/",sample_type[s],"_SurfaceGenie_marker_gene_table_subset_arranged.rds")) # #Level2_
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
        group <- paste(datasets[d],sample_type[s],sep = ":")
        markers2 <- as.data.frame(cbind (markers=markers,group = rep(group,length(markers))))
        markers3 <- as.data.frame(cbind (markers=neg_markers,group = rep(group,length(neg_markers))))
        markers_chose <- rbind(markers_chose,markers2)
        neg_markers_chose <- rbind(neg_markers_chose,markers3)
      }
      
      
    }
    
    markers_chose$value <- rep(1,dim(markers_chose)[1])
    neg_markers_chose$value <- rep(-1,dim(neg_markers_chose)[1])
    heatmap_pos <- as.data.frame(tidyr::pivot_wider(markers_chose, names_from = group, values_from = value))
    heatmap_neg <- as.data.frame(tidyr::pivot_wider(neg_markers_chose, names_from = group, values_from = value))
    
    heatmap <- rbind(heatmap_pos,heatmap_neg)
    heatmap[is.na(heatmap)] <- 0
    heatmap <- as.data.frame(heatmap)
    col_order <- c("markers","travaglini_2020:general_normal","bischoff_2021:adjacent_normal","kim_2020:adjacent_normal","lambrechts_2018:adjacent_normal", "song_2019:adjacent_normal",
                   "bischoff_2021:luad","lambrechts_2018:luad","wu_2021:luad",
                   "lambrechts_2018:lusc" ,"wu_2021:lusc",
                   "kim_2020:nsclc_primary","lambrechts_2018:nsclc_primary","wu_2021:nsclc_primary","zilionis_2019:nsclc_primary" ,"song_2019:nsclc_primary",
                   "kim_2020:nsclc_meta",
                   "kim_2020:pleural_fluids","kim_2020:LymphNode_normal", "kim_2020:LymphNode_meta",
                   "zilionis_2019:blood")
    heatmap <-aggregate(heatmap[,-1], by=list(markers=heatmap$markers), FUN=sum)
    
    heatmap <-heatmap[,col_order]
    row.names.remove <- c("No markers selected")
    heatmap <-heatmap[!(heatmap$markers %in% row.names.remove),]
    club_colname <- paste(col_order,"Club")
    AT2_colname <- paste(col_order,"AT2")
    
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
  annotation <- as.data.frame(stringr::str_split_fixed(df$x, ":", 2))
  annotation2 <- as.data.frame(stringr::str_split_fixed(annotation$V2, " ", 2))
  annotation <- cbind(annotation,annotation2)
  annotation <- annotation[,c(1,3,4)]
  row.names(annotation) <- colnames(Sheatmap)
  colnames(annotation) <- c("Cohort","Sample Type","Cell Type")
  Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
           "#661100", "#0072B2", "#D55E00")
  names(Var1) = unique(sort(annotation$Cohort))
  Var2 = c("lightblue", "darkred","navy","pink2","purple",
           "gold","orange","greenyellow","darkgreen","royalblue")
  names(Var2) = unique(sort(annotation$"Sample Type"))
  Var3 <- c("moccasin", "firebrick")
  names(Var3) = c("Club","AT2")
  
  ann_colors = list(Cohort = Var1, "Sample Type" = Var2,"Cell Type"= Var3)
  dim(Sheatmap) #189 40
  Sheatmap_copy <- Sheatmap
  Sheatmap_copy$rowSum <- rowSums(Sheatmap_copy)
  
  Sheatmap1 <- Sheatmap_copy[Sheatmap_copy$rowSum>0,1:40]
  Sheatmap2 <- Sheatmap_copy[Sheatmap_copy$rowSum<0,1:40]
  
  pdf(paste0(outdir,"/20220202_Heatmap_SurfaceGenie_marker_chose_clusterRow_AT2_Club_compare.pdf"), width=12, height=30)
  
  a <- pheatmap::pheatmap(Sheatmap,
                          main = "SurfaceGenie Markers, public datasets",
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
  
  pdf(paste0(outdir,"/20220202_Heatmap_SurfaceGenie_marker_chose_clusterRow_AT2_Club_separate.pdf"), width=12, height=25)
  
  a <- pheatmap::pheatmap(Sheatmap1,
                          main = "SurfaceGenie Markers -1, public datasets",
                          color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                          annotation = annotation,
                          legend_breaks= c(1,0,-1),
                          legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                          annotation_colors = ann_colors, 
                          cluster_rows = T,
                          cellheight=10,
                          cluster_cols = F,
                          show_colnames= F)
  
  print(a)
  b <- pheatmap::pheatmap(Sheatmap2,
                          main = "SurfaceGenie Markers -2, public datasets",
                          color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                          annotation = annotation,
                          legend_breaks= c(1,0,-1),
                          legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                          annotation_colors = ann_colors, 
                          cluster_rows = T,
                          cellheight=10,
                          cluster_cols = F,
                          show_colnames= F)
  print(b)
  dev.off()

## SurfaceGenie Vs GO #######
intersect2 <- intersect(row.names(heatmap),row.names(Sheatmap))
new_heatmap <- Sheatmap[!(row.names(Sheatmap)%in%intersect2 ),]
pdf(paste0(outdir,"/20220202_Heatmap_NewS_marker_chose_clusterRow_AT2_Club_compare.pdf"), width=12, height=25)

a <- pheatmap::pheatmap(new_heatmap,
                        main = "New SurfaceGenie Markers, public datasets",
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
c <- 1L
for(c in 1:length(cell_in_interesed)){
  cell_type <- cell_in_interesed[c]
  print(cell_type)
  
  
  markers_chose <- c()
  neg_markers_chose <- c()
  
  #s <- 5L 
  for( s in 1:length(sample_type)){
    disease <- sample_type[s]
    print(disease)
    datasets <- list_sample_cohort[[disease]]
    
    #d <- 1L
    for(d in 1:length(datasets)){
      print(datasets[d])
      ## read files ####
      table <- readRDS(paste0(outdir,"/",datasets[d],"/Level2_",sample_type[s],"_marker_gene_table_subset_arranged.rds")) #_SurfaceGenie
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
      group <- paste(datasets[d],sample_type[s],sep = ":")
      markers2 <- as.data.frame(cbind (markers=markers,group = rep(group,length(markers))))
      markers3 <- as.data.frame(cbind (markers=neg_markers,group = rep(group,length(neg_markers))))
      markers_chose <- rbind(markers_chose,markers2)
      neg_markers_chose <- rbind(neg_markers_chose,markers3)
    }
    
    
  }
  
  markers_chose$value <- rep(1,dim(markers_chose)[1])
  neg_markers_chose$value <- rep(-1,dim(neg_markers_chose)[1])
  heatmap_pos <- as.data.frame(tidyr::pivot_wider(markers_chose, names_from = group, values_from = value))
  heatmap_neg <- as.data.frame(tidyr::pivot_wider(neg_markers_chose, names_from = group, values_from = value))
  
  heatmap <- rbind(heatmap_pos,heatmap_neg)
  
  heatmap[is.na(heatmap)] <- 0
  heatmap <- as.data.frame(heatmap)
  col_order <- c("markers","travaglini_2020:general_normal","bischoff_2021:adjacent_normal","kim_2020:adjacent_normal","lambrechts_2018:adjacent_normal", "song_2019:adjacent_normal",
                 "bischoff_2021:luad","lambrechts_2018:luad","wu_2021:luad",
                 "lambrechts_2018:lusc" ,"wu_2021:lusc",
                 "kim_2020:nsclc_primary","lambrechts_2018:nsclc_primary","wu_2021:nsclc_primary","zilionis_2019:nsclc_primary" ,"song_2019:nsclc_primary",
                 "kim_2020:nsclc_meta",
                 "kim_2020:pleural_fluids","kim_2020:LymphNode_normal", "kim_2020:LymphNode_meta",
                 "zilionis_2019:blood")
  heatmap <-aggregate(heatmap[,-1], by=list(markers=heatmap$markers), FUN=sum)
  heatmap <-heatmap[,col_order]
  row.names.remove <- c("No markers selected")
  heatmap <-heatmap[!(heatmap$markers %in% row.names.remove),]
  club_colname <- paste(col_order,"Club")
  AT2_colname <- paste(col_order,"AT2")
  
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
annotation <- as.data.frame(stringr::str_split_fixed(df$x, ":", 2))
annotation2 <- as.data.frame(stringr::str_split_fixed(annotation$V2, " ", 2))
annotation <- cbind(annotation,annotation2)
annotation <- annotation[,c(1,3,4)]
row.names(annotation) <- colnames(heatmap)
colnames(annotation) <- c("Cohort","Sample Type","Cell Type")
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
         "#661100", "#0072B2", "#D55E00")
names(Var1) = unique(sort(annotation$Cohort))
Var2 = c("lightblue", "darkred","navy","pink2","purple",
         "gold","orange","greenyellow","darkgreen","royalblue")
names(Var2) = unique(sort(annotation$"Sample Type"))
Var3 <- c("moccasin", "firebrick")
names(Var3) = c("Club","AT2")

ann_colors = list(Cohort = Var1, "Sample Type" = Var2,"Cell Type"= Var3)

pdf(paste0(outdir,"/20220202_Level2_Heatmap_GOmarker_chose_clusterRow_AT2_Club_compare.pdf"), width=12, height=20)

a <- pheatmap::pheatmap(heatmap,
                        main = "GO Markers, public datasets (Level 2)",
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
  
  s <- 3L 
  for( s in 1:length(sample_type)){
    disease <- sample_type[s]
    print(disease)
    datasets <- list_sample_cohort[[disease]]
    
    #d <- 1L
    for(d in 1:length(datasets)){
      print(datasets[d])
      ## read files ####
      table <- readRDS(paste0(outdir,"/",datasets[d],"/Level2_",sample_type[s],"_SurfaceGenie_marker_gene_table_subset_arranged.rds")) # #Level2_
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
      group <- paste(datasets[d],sample_type[s],sep = ":")
      markers2 <- as.data.frame(cbind (markers=markers,group = rep(group,length(markers))))
      markers3 <- as.data.frame(cbind (markers=neg_markers,group = rep(group,length(neg_markers))))
      markers_chose <- rbind(markers_chose,markers2)
      neg_markers_chose <- rbind(neg_markers_chose,markers3)
    }
    
    
  }
  
  markers_chose$value <- rep(1,dim(markers_chose)[1])
  neg_markers_chose$value <- rep(-1,dim(neg_markers_chose)[1])
  heatmap_pos <- as.data.frame(tidyr::pivot_wider(markers_chose, names_from = group, values_from = value))
  heatmap_neg <- as.data.frame(tidyr::pivot_wider(neg_markers_chose, names_from = group, values_from = value))
  
  heatmap <- rbind(heatmap_pos,heatmap_neg)
  heatmap[is.na(heatmap)] <- 0
  heatmap <- as.data.frame(heatmap)
  col_order <- c("markers","travaglini_2020:general_normal","bischoff_2021:adjacent_normal","kim_2020:adjacent_normal","lambrechts_2018:adjacent_normal", "song_2019:adjacent_normal",
                 "bischoff_2021:luad","lambrechts_2018:luad","wu_2021:luad",
                 "lambrechts_2018:lusc" ,"wu_2021:lusc",
                 "kim_2020:nsclc_primary","lambrechts_2018:nsclc_primary","wu_2021:nsclc_primary","zilionis_2019:nsclc_primary" ,"song_2019:nsclc_primary",
                 "kim_2020:nsclc_meta",
                 "kim_2020:pleural_fluids","kim_2020:LymphNode_normal", "kim_2020:LymphNode_meta",
                 "zilionis_2019:blood")
  heatmap <-aggregate(heatmap[,-1], by=list(markers=heatmap$markers), FUN=sum)
  
  heatmap <-heatmap[,col_order]
  row.names.remove <- c("No markers selected")
  heatmap <-heatmap[!(heatmap$markers %in% row.names.remove),]
  club_colname <- paste(col_order,"Club")
  AT2_colname <- paste(col_order,"AT2")
  
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
annotation <- as.data.frame(stringr::str_split_fixed(df$x, ":", 2))
annotation2 <- as.data.frame(stringr::str_split_fixed(annotation$V2, " ", 2))
annotation <- cbind(annotation,annotation2)
annotation <- annotation[,c(1,3,4)]
row.names(annotation) <- colnames(Sheatmap)
colnames(annotation) <- c("Cohort","Sample Type","Cell Type")
Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
         "#661100", "#0072B2", "#D55E00")
names(Var1) = unique(sort(annotation$Cohort))
Var2 = c("lightblue", "darkred","navy","pink2","purple",
         "gold","orange","greenyellow","darkgreen","royalblue")
names(Var2) = unique(sort(annotation$"Sample Type"))
Var3 <- c("moccasin", "firebrick")
names(Var3) = c("Club","AT2")

ann_colors = list(Cohort = Var1, "Sample Type" = Var2,"Cell Type"= Var3)
dim(Sheatmap) #151
Sheatmap_copy <- Sheatmap
Sheatmap_copy$rowSum <- rowSums(Sheatmap_copy)

Sheatmap1 <- Sheatmap_copy[Sheatmap_copy$rowSum>0,1:40]
Sheatmap2 <- Sheatmap_copy[Sheatmap_copy$rowSum<0,1:40]

pdf(paste0(outdir,"/20220202_Level2_Heatmap_SurfaceGenie_marker_chose_clusterRow_AT2_Club_compare.pdf"), width=12, height=30)

a <- pheatmap::pheatmap(Sheatmap,
                        main = "SurfaceGenie Markers, public datasets (Level 2)",
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

pdf(paste0(outdir,"/20220202_Level2_Heatmap_SurfaceGenie_marker_chose_clusterRow_AT2_Club_separate.pdf"), width=12, height=14)

a <- pheatmap::pheatmap(Sheatmap1,
                        main = "SurfaceGenie Markers -1, public datasets (Level 2)",
                        color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                        annotation = annotation,
                        legend_breaks= c(1,0,-1),
                        legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                        annotation_colors = ann_colors, 
                        cluster_rows = T,
                        cellheight=10,
                        cluster_cols = F,
                        show_colnames= F)

print(a)
b <- pheatmap::pheatmap(Sheatmap2,
                        main = "SurfaceGenie Markers -2, public datasets (Level 2)",
                        color = colorRampPalette(c("steelblue","white", "firebrick"))(3),
                        annotation = annotation,
                        legend_breaks= c(1,0,-1),
                        legend_labels= c("Chosen pos.", "Not Chosen","Chosen neg."),
                        annotation_colors = ann_colors, 
                        cluster_rows = T,
                        cellheight=10,
                        cluster_cols = F,
                        show_colnames= F)
print(b)
dev.off()

## SurfaceGenie Vs GO #######
intersect2 <- intersect(row.names(heatmap),row.names(Sheatmap))
new_heatmap <- Sheatmap[!(row.names(Sheatmap)%in%intersect2 ),]
pdf(paste0(outdir,"/20220202_Level2_Heatmap_NewS_marker_chose_clusterRow_AT2_Club_compare.pdf"), width=12, height=22)

a <- pheatmap::pheatmap(new_heatmap,
                        main = "New SurfaceGenie Markers, public datasets (Level 2)",
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
features <- rownames(Sheatmap) #using non-stratified #189

heatmap_data <- data.frame("gene"=features)
empty <- data.frame("gene"=features)

#f <- 1L
for( f in 1:length(filename)){
  print(filename[f])
  ob <- readRDS(paste0(rawdir,"/",filename[f],".rds"))
  
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
      RotatedAxis()+  labs(title=paste0(filename[f],": ",split[k]))+
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
      colnames(data_cell) <- c(paste0(cell_type[c],":",filename[f],":",split[k]),"gene")
      data2 <- merge(data2,data_cell,by="gene",all=T) 
    }
    
    
    
    heatmap_data <- merge(heatmap_data,data2,by="gene",all=T)
    rm(ob_split_PTPRC0)
    rm(ob2) 
    rm(ob22)
  }  
  
}
#heatmap_data2 <- heatmap_data[, colSums(is.na(heatmap_data)) != nrow(heatmap_data)]
#heatmap_data2 <-heatmap_data %>% select_if(~sum(!is.na(.)) > 0)
heatmap_data[is.na(heatmap_data)] <- 0
heatmap_data <- as.data.frame(heatmap_data)
row.names(heatmap_data) <- heatmap_data$"gene"
heatmap_data$"gene"<- NULL
saveRDS(heatmap_data,paste0(outdir,"/20220202_Heatmap_SurfaceGenie_expression_data_level2.rds"))

## can start here!!!!!!!#############
heatmap_data <- readRDS(paste0(outdir,"/20220202_Heatmap_SurfaceGenie_expression_data_level2.rds"))
df <- data.frame(x = colnames(heatmap_data))
annotation <- as.data.frame(stringr::str_split_fixed(df$x, ":", 2))
annotation2 <- as.data.frame(stringr::str_split_fixed(annotation$V2, ":", 2))
annotation <- cbind(annotation,annotation2)
annotation <- annotation[,c(1,3,4)]
row.names(annotation) <- colnames(heatmap_data)
colnames(annotation) <- c("Cell Type","Cohort","Sample Type")
annotation$`Cell Type` <- gsub("Signaling Alveolar Epithelial Type 2","Alveolar Epithelial Type 2",annotation$`Cell Type`)

Var1 = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
         "#661100", "#0072B2", "#D55E00")
names(Var1) = unique(sort(annotation$Cohort))
Var2 = c("lightblue", "darkred","navy","pink2","purple",
         "gold","orange","greenyellow","darkgreen","royalblue")
names(Var2) = unique(sort(annotation$"Sample Type"))
Var3 <- c("hotpink", "firebrick","green2","gold","moccasin","darkgreen","black",
          "steelblue1","gray","navy","limegreen","seagreen2","orange","purple1")
names(Var3) = unique(sort(annotation$"Cell Type"))


ann_colors = list(Cohort = Var1, "Sample Type" = Var2,"Cell Type"= Var3)


### non- pearson ##########
cell_order <- c("Alveolar Epithelial Type 2" ,
                "Club","Alveolar Epithelial Type 1",
                "Mucous" ,"Goblet","Basal","Differentiating Basal","Proliferating Basal",
                "Proximal Basal", 
                "Ciliated" ,"Proximal Ciliated","Ionocyte","Neuroendocrine"  ,"Serous")
heatmap_sub <- heatmap_data %>% dplyr:: select(starts_with(cell_order))


pdf(paste0(outdir,"/20220202_Heatmap_SurfaceGenie_experssion_clusterRow_euclidean.pdf"), width=16, height=28)

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

pdf(paste0(outdir,"/20220202_Heatmap_SurfaceGenie_experssion_clusterRow_Col_euclidean.pdf"), width=16, height=28)

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

