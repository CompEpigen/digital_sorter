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
outdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/lung_mapped_cellxgene_fixed_addsplits"


# Load functions
GetfileNames <- function(fileDir, pattern = ".rds"){
  filenames <- list.files(fileDir, pattern= pattern, full.names=TRUE)
  filename <- sub(paste0(".*",fileDir,"/"), "", filenames)
  filename <- sub(pattern, "", filename)  
  return(filename)
}


# Read files ####
cell.surface.marker <-readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/digital_sorter/R/cell.surface.marker.rds")
features <- c("SCGB1A1","UPK3A","SCGB3A2","MET","MGP","IVL","UPK2","CTSW", #Possible club-related genes
              "CD74","TFPI","SDC4", #In both club and AT2 cells
              "SFTPA1","SFTPA2","SFTPB","SFTPC","SFTPD","FABP5","CBR2","LGI3","ABCA3","LAMP3", #AT2 
              "LCN2","IL33" ,"IFI27L2A", #activated AT2 before transitioning into KRT8+ intermediate cells #COPD
              "KRT8","PLAUR","SPRR1A","CDKN1A","CLDN4","AREG","HBEGF","ITGB6","EDN1", #KRT8+ cells
              "RTKN2","PDPN","HOPX","AGER","IGFBP2" #AT1
)

filename <- GetfileNames(rawdir, pattern = ".rds")

f <- 6L
for( f in 1:length(filename)){
filename[f]
ob <- readRDS(paste0(rawdir,"/",filename[f],".rds"))
  
#dir.create(paste0(outdir,"/",filename))
split <- as.character(unique(sort(ob@meta.data$disease)))
ob_split <- SplitObject(ob, split.by = "disease")
length(split)

k <- 3L
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

pdf(paste0(outdir,"/",filename[f],"/Dotplot_",filename[f],"_",split[k],"_CD45-EPCAM+_GenesInterested2.pdf"), width=20, height=6)
DotPlot(ob_split_PTPRC0, features = features, group.by = "annotation.l2") + 
  RotatedAxis()+  labs(title=paste0(filename[f],": ",split[k]))+
  theme(plot.title = element_text(hjust = 0.5)) 
dev.off()
rm(ob_split_PTPRC0)
rm(ob2) 
rm(ob22)
}  

}   
    