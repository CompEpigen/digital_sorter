library(dplyr)
library("viridis") 
library(ggplot2)
library(UpSetR)
library(Seurat)
library(SeuratObject)
# Set working directory ####
wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)
rawdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/markerlist_for_heatmap"
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

epithelial_cells <-unique(sort(ob_split2@meta.data[["annotation.l2"]]))


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
 

 # show only gene_interest
 source("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/markerlist_for_heatmap/util.R")
 
 GeneFilesnames <- GetfileNames(rawdir)
 i<-3L
 for(i in seq(length(GeneFilesnames))){
   geneFiles <- read.csv(paste0(rawdir,"/",GeneFilesnames[i],".csv"))
   genes <-geneFiles[,1,drop=T]
   Sheatmap2 <- Sheatmap%>% filter(rownames(Sheatmap)%in% genes)
   Sheatmap2 <- Sheatmap2[genes,]
   Sheatmap2 <- na.omit(Sheatmap2)
   
   pdf(paste0(outdir,"/level2_marker_selection_table/Heatmap_",GeneFilesnames[i],"_epithelial.pdf"), width=14, height=8)
   
   a <- pheatmap::pheatmap(Sheatmap2,
                           main = paste(GeneFilesnames[i],", our own dataset with tumor annotation (Level2)"),
                           
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
##heatmap expresion#####
 features <- genes ###dependent to gene files
 
 heatmap_data <- data.frame("gene"=features)
 empty <- data.frame("gene"=features)
 
 
 ob22 <- SplitObject(ob_split2, split.by = "disease")
 
 
 k <- 1L
 for(k in 1:length(sample_type)){
   print(sample_type[k])      
   ob_split_PTPRC0 = ob22[[sample_type[k]]]
   
   
   a <- DotPlot(ob_split_PTPRC0, features = features, group.by = "annotation.l2") + 
     RotatedAxis()+  
     theme(plot.title = element_text(hjust = 0.5)) 
   data <- a[["data"]]
   data <- data[,c(1,3,4)] #avg.exp
   cell_type <- unique(data$id)
   data2 <- empty
   #c <- 2L
   for(c in seq(length(cell_type))){
     print(cell_type[c])
     data_cell <- data %>% filter(id == cell_type[c])
     data_cell$id <- NULL
     colnames(data_cell) <- c(paste0(cell_type[c],":",sample_type[k]),"gene")
     data2 <- merge(data2,data_cell,by="gene",all=T) 
   }
   
   
   
   heatmap_data <- merge(heatmap_data,data2,by="gene",all=T)
   rm(ob_split_PTPRC0)
   
 }    
 
 heatmap_data[is.na(heatmap_data)] <- 0
 heatmap_data <- as.data.frame(heatmap_data)
 row.names(heatmap_data) <- heatmap_data$"gene"
 heatmap_data$"gene"<- NULL
 
 
 df <- data.frame(x = colnames(heatmap_data))
 annotation <- as.data.frame(stringr::str_split_fixed(df$x, ":", 2))
 row.names(annotation) <- colnames(heatmap_data)
 colnames(annotation) <- c("Cell Type","Sample Type")
 annotation$`Cell Type` <- gsub("Signaling Alveolar Epithelial Type 2","Alveolar Epithelial Type 2",annotation$`Cell Type`)
 
 
 Var3 <- c("hotpink", "firebrick","green2",
           "blue2","royalblue4","steelblue4","deepskyblue",
           "grey","gold","moccasin","darkgreen","black",
           "lightblue","orangered","navy","limegreen","turquoise","orange","purple1")
 names(Var3) = unique(sort(annotation$"Cell Type"))
 names(Var1) = unique(sort(annotation$"Sample Type"))[c(2,1,4,3)]
 
 ann_colors = list("Sample Type" = Var1,"Cell Type"= Var3)
 
 
 cell_order <- c("C9_Unknown", "C1_Tumor","C4_Tumor","C5_Tumor","C15_Tumor",
                 "Signaling Alveolar Epithelial Type 2","Alveolar Epithelial Type 2" ,
                 "Club","Alveolar Epithelial Type 1",
                 "Mucous" ,"Goblet","Basal","Differentiating Basal","Proliferating Basal",
                 "Proximal Basal", 
                 "Ciliated" ,"Proximal Ciliated","Ionocyte","Neuroendocrine"  ,"Serous")
 heatmap_sub <- heatmap_data %>% dplyr:: select(starts_with(cell_order))
 #heatmap_sub <- na.omit(heatmap_sub)
 
 heatmap_sub2 <-heatmap_sub[rowSums(heatmap_sub == 0) < 63, ]
 
 pdf(paste0(outdir,"/Heatmap_expression",GeneFilesnames[i],"_epithelial.pdf"), width=13, height=31)

 a <- pheatmap::pheatmap(heatmap_sub2,
                         color = colorRampPalette(c("blue4","white", "red2"))(50),
                         main = paste(GeneFilesnames[i]),
                         
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
 
 
 a <- pheatmap::pheatmap(heatmap_sub2,
                         color = colorRampPalette(c("blue4","white", "red2"))(50),
                         main = paste(GeneFilesnames[i]),
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
 }
 
 