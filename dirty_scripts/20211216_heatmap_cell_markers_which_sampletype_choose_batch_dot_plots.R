library(dplyr)
library("viridis") 
library(ggplot2)
library(UpSetR)
# Set working directory ####
wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)
rawdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/lung_mapped_cellxgene_fixed_addsplits"
outdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/lung_mapped_cellxgene_fixed_addsplits"

## Create list ####
if(F){ 
  sample_type <- c("general_normal","adjacent_normal","blood","luad","lusc",
                   "nsclc_primary","nsclc_meta","LymphNode_normal","LymphNode_meta","pleural_fluids")
  general_normal  <- c("travaglini_2020")
  adjacent_normal <- c("song_2019","lambrechts_2018","bischoff_2021","kim_2020" )
  blood           <- c("zilionis_2019")
  luad            <- c("lambrechts_2018","wu_2021","bischoff_2021")
  lusc            <- c("lambrechts_2018","wu_2021")
  nsclc_primary   <- c("song_2019","lambrechts_2018","zilionis_2019","wu_2021","kim_2020" )
  nsclc_meta      <- c("kim_2020")
  LymphNode_normal<- c("kim_2020")
  LymphNode_meta  <- c("kim_2020")
  pleural_fluids  <- c("kim_2020")
  
  list_sample_cohort <- list(general_normal, adjacent_normal, blood, luad, lusc, nsclc_primary,
                             nsclc_meta, LymphNode_normal, LymphNode_meta , pleural_fluids)
  names(list_sample_cohort) <- sample_type
  saveRDS(list_sample_cohort, paste0(outdir,"/list_sample_cohort.rds"))
}


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


c <- 1L
for(c in 1:length(cell_in_interesed)){
  cell_type <- cell_in_interesed[c]
  print(cell_type)
  
  
  markers_chose <- c()
  
  #s <- 10L #1-10
  for( s in 1:length(sample_type)){
    disease <- sample_type[s]
    print(disease)
    datasets <- list_sample_cohort[[disease]]
    
    
    
    #d <- 1L
    for(d in 1:length(datasets)){
      print(datasets[d])
      ## read files ####
      table <- readRDS(paste0(outdir,"/",datasets[d],"/",sample_type[s],"_marker_gene_table_subset_arranged20_exp.rds"))
      table <- table %>% filter(table$annotation.l2 == cell_type)
      markers <- table$gene
      markers <- unique(sort(markers))
      if(length(markers)==0){
        markers <- "No markers selected"
      }
      group <- paste(datasets[d],sample_type[s],sep = ":")
      markers2 <- as.data.frame(cbind (markers=markers,group = rep(group,length(markers))))
      markers_chose <- rbind(markers_chose,markers2)
    }
    
    
  }

  markers_chose$value <- rep(1,dim(markers_chose)[1])
  heatmap <- as.data.frame(tidyr::pivot_wider(markers_chose, names_from = group, values_from = value))
  col_order <- c("markers","travaglini_2020:general_normal","bischoff_2021:adjacent_normal","kim_2020:adjacent_normal","lambrechts_2018:adjacent_normal", "song_2019:adjacent_normal",
                 "bischoff_2021:luad","lambrechts_2018:luad","wu_2021:luad",
                 "lambrechts_2018:lusc" ,"wu_2021:lusc",
                 "kim_2020:nsclc_primary","lambrechts_2018:nsclc_primary","wu_2021:nsclc_primary","zilionis_2019:nsclc_primary" ,"song_2019:nsclc_primary",
                 "kim_2020:nsclc_meta",
                 "kim_2020:pleural_fluids","kim_2020:LymphNode_normal", "kim_2020:LymphNode_meta",
                 "zilionis_2019:blood")
  heatmap[is.na(heatmap)] <- 0
  heatmap <- as.data.frame(heatmap)
  
  #row.names(heatmap) <- heatmap$markers
  #heatmap$markers <- NULL
  heatmap <-heatmap[,col_order]
  row.names.remove <- c("No markers selected")
  heatmap <-heatmap[!(heatmap$markers %in% row.names.remove),]
  club_colname <- paste(col_order,"Club")
  AT2_colname <- paste(col_order,"AT2")
  # Club
  heatmap_club <- heatmap
  colnames(heatmap_club) <- club_colname
  
  # AT2
  heatmap_AT2 <- heatmap
  colnames(heatmap_AT2) <- AT2_colname
  
  heatmap <- merge(heatmap_club,heatmap_AT2,
                   by.x = "markers Club",
                   by.y = "markers AT2",all =T)
  
  row.names(heatmap) <- heatmap$"markers Club"
  #######################################
  #get AT2_exclude and club_exclude from the upset list!!!!
  c <- 2L #change here
  
    cell_type <- cell_in_interesed[c]
    print(cell_type)
    
    upset_list <- list()
    
    #s <- 4L #1-10
    for( s in 1:length(sample_type)){
      disease <- sample_type[s]
      print(disease)
      datasets <- list_sample_cohort[[disease]]
      
      
      union_markers <- c()
      #d <- 2L
      for(d in 1:length(datasets)){
        print(datasets[d])
        ## read files ####
        table <- readRDS(paste0(outdir,"/",datasets[d],"/",sample_type[s],"_marker_gene_table_subset_arranged20_exp.rds"))
        table <- table %>% filter(table$annotation.l2 == cell_type)
        markers <- table$gene
        markers <- unique(sort(markers))
        union_markers <- append(union_markers,markers)
      }
      union_markers <-unique(sort(union_markers))
      list <- list(union_markers)
      names(list)<- disease
      upset_list <- append(upset_list,list)
      
    }
    count <- table(unlist(upset_list))
    for_barplot <- data.frame(cbind(count))
    for_barplot$gene <- row.names(for_barplot)
    
  if(c == 1L){ #when c <- 1L
    club_exclude <- for_barplot %>% filter(for_barplot$count ==1)
    club_exclude <- club_exclude$gene
  }
  if(c == 2L){ #when c <- 2L
    AT2_exclude <- for_barplot %>% filter(for_barplot$count ==1)
    AT2_exclude <- AT2_exclude$gene
  }
    
    ###########################################
  heatmap <-heatmap %>% filter(!(heatmap$"markers Club"%in%c(AT2_exclude,club_exclude)))
  
  heatmap$"markers Club" <- NULL
  heatmap[is.na(heatmap)] <- 0
 
  df <- data.frame(x = colnames(heatmap))
  annotation <- as.data.frame(stringr::str_split_fixed(df$x, ":", 2))
  annotation2 <- as.data.frame(stringr::str_split_fixed(annotation$V2, " ", 2))
  annotation <- cbind(annotation,annotation2)
  annotation <- annotation[,c(1,3,4)]
  row.names(annotation) <- colnames(heatmap)
  colnames(annotation) <- c("Cohort","Sample Type","Cell Type")
  Var1 = rainbow(7)
  names(Var1) = unique(sort(annotation$Cohort))
  Var2 = c("lightblue", "darkred","navy","pink2","purple",
           "gold","orange","greenyellow","darkgreen","royalblue")
  names(Var2) = unique(sort(annotation$"Sample Type"))
  Var3 <- c("moccasin", "firebrick")
  names(Var3) = c("Club","AT2")
  
  
  ann_colors = list(Cohort = Var1, "Sample Type" = Var2,"Cell Type"= Var3)
  
  pdf(paste0(outdir,"/Heatmap_marker_chose_clusterRow_",cell_type,"_compare.pdf"), width=10, height=10)
  
  a <- pheatmap::pheatmap(heatmap,
                          color = colorRampPalette(c("white", "red"))(50),
                     main = "Markers chosen by at least 2 cohorts",
                     annotation = annotation,
                     annotation_colors = ann_colors, 
                     cluster_rows = T,
                     cluster_cols = F)
  print(a)
  
  dev.off()
}


# drawing dotplots for CD45- EPCAM+ cells ####
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)


# Set working directory ####
wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)
rawdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_each_dataset"
outdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/lung_mapped_cellxgene_fixed_addsplits/1216Dotplots"


# Load functions
GetfileNames <- function(fileDir, pattern = ".rds"){
  filenames <- list.files(fileDir, pattern= pattern, full.names=TRUE)
  filename <- sub(paste0(".*",fileDir,"/"), "", filenames)
  filename <- sub(pattern, "", filename)  
  return(filename)
}


# Read files ####
cell.surface.marker <-readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/digital_sorter/R/cell.surface.marker.rds")
#according to heatmap
features <- c("B2M", "BCAM" ,"CEACAM6","EPCAM","KRT18", "SDC4", #shared markers
              "AQP4","CLU","FOLR1" ,"LMO7", #shared markers (not shared in many sample types)
              "MMP7" ,"ADGRF5" ,"C3","CTSB","LGALS1","SDC1","EMP2" ,"TFPI" , #shared (not shared in normal)
              "HHIP" , "PLA2G1B" ,"WIF1","CSF3R", "SERPINF2",  #AT2
              "AGER" ,"ANXA2","HPN"  #club  
              )

filename <- GetfileNames(rawdir, pattern = ".rds")

#f <- 4L
for( f in 1:length(filename)){
  print(filename[f])
  ob <- readRDS(paste0(rawdir,"/",filename[f],".rds"))
  
  #dir.create(paste0(outdir,"/",filename))
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
    
    # plot dot plot
    ob_split_PTPRC0@meta.data[["annotation.l2"]] <- factor(ob_split_PTPRC0@meta.data[["annotation.l2"]], 
                                                           levels=c("Signaling Alveolar Epithelial Type 2", "Alveolar Epithelial Type 2",
                                                                    "Club" ,"Alveolar Epithelial Type 1",
                                                                    "Mucous" , "Goblet" ,"Basal", "Differentiating Basal", "Proximal Basal" ,"Proliferating Basal", 
                                                                    "Ionocyte" ,"Neuroendocrine" ,
                                                                    "Ciliated", "Proximal Ciliated" ,"Serous"))
    
    pdf(paste0(outdir,"/Dotplot_",filename[f],"_",split[k],"_CD45-EPCAM+_Genes2.pdf"), width=20, height=6)
    a <- DotPlot(ob_split_PTPRC0, features = features, group.by = "annotation.l2") + 
      RotatedAxis()+  labs(title=paste0(filename[f],": ",split[k]))+
      theme(plot.title = element_text(hjust = 0.5)) 
    print(a)
    dev.off()
    rm(ob_split_PTPRC0)
    rm(ob2) 
    rm(ob22)
  }  
  
}   



