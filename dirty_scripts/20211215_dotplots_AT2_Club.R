
# drawing dotplots for CD45- EPCAM+ cells ####
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)


# Set working directory ####
wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/lung_mapped_cellxgene_fixed_addsplits/1_upsetplot_cellType_Wenn"
setwd(wd)
rawdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_each_dataset"
outdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/lung_mapped_cellxgene_fixed_addsplits/1215Dotplots"


# Load functions
GetfileNames <- function(fileDir, pattern = ".rds"){
  filenames <- list.files(fileDir, pattern= pattern, full.names=TRUE)
  filename <- sub(paste0(".*",fileDir,"/"), "", filenames)
  filename <- sub(pattern, "", filename)  
  return(filename)
}

# Read files ####
wenn_list <- readRDS("wenn_list.rds")
ItemsList <- gplots::venn(wenn_list, show.plot = FALSE)
filename <- GetfileNames(rawdir, pattern = ".rds")
# features
Club <- attributes(ItemsList)$intersections[["Club"]]
AT2 <- attributes(ItemsList)$intersections[["Alveolar Epithelial Type 2"]]
Share <-attributes(ItemsList)$intersections[["Club:Alveolar Epithelial Type 2"]]


#f <- 7L
for( f in 1:length(filename)){
  print(filename[f])
  ob <- readRDS(paste0(rawdir,"/",filename[f],".rds"))
  
  #dir.create(paste0(outdir,"/",filename))
  split <- as.character(unique(sort(ob@meta.data$disease)))
  ob_split <- SplitObject(ob, split.by = "disease")
  length(split)
  
  k <- 2L
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
    
    pdf(paste0(outdir,"/Dotplot_",filename[f],"_",split[k],"_CD45-EPCAM+_GenesInterested_AT2.pdf"), width=20, height=6)
    a <- DotPlot(ob_split_PTPRC0, features = AT2, group.by = "annotation.l2") + 
      RotatedAxis()+  labs(title=paste0(filename[f],": ",split[k]," AT2"))+
      theme(plot.title = element_text(hjust = 0.5)) 
    print(a)
    dev.off()
    pdf(paste0(outdir,"/Dotplot_",filename[f],"_",split[k],"_CD45-EPCAM+_GenesInterested_Club.pdf"), width=20, height=6)
    a <- DotPlot(ob_split_PTPRC0, features = Club, group.by = "annotation.l2") + 
      RotatedAxis()+  labs(title=paste0(filename[f],": ",split[k]," Club"))+
      theme(plot.title = element_text(hjust = 0.5)) 
    print(a)
    dev.off()
    pdf(paste0(outdir,"/Dotplot_",filename[f],"_",split[k],"_CD45-EPCAM+_GenesInterested_Share.pdf"), width=20, height=6)
    a <- DotPlot(ob_split_PTPRC0, features = Share, group.by = "annotation.l2") + 
      RotatedAxis()+  labs(title=paste0(filename[f],": ",split[k]," Share"))+
      theme(plot.title = element_text(hjust = 0.5)) 
    print(a)
    dev.off()
    rm(ob_split_PTPRC0)
    rm(ob2) 
    rm(ob22)
  }  
  
}   
