library(dplyr)
library("viridis") 
library(ggplot2)
library(UpSetR)
# Set working directory ####
wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)
rawdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/patient_hlcma_addsplits"
outdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/patient_hlcma_addsplits"

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

sample_type <- c("distal_normal","luadc","tan")

cell_in_interesed <- c("Club","Alveolar Epithelial Type 2")
wenn_list <- list()

c <- 1L
for(c in 1:length(cell_in_interesed)){
  cell_type <- cell_in_interesed[c]
  print(cell_type)
  
  upset_list <- list()
  
  #s <- 4L #1-10
  for( s in 1:length(sample_type)){
    disease <- sample_type[s]
    print(disease)
    
    table <- readRDS(paste0(outdir,"/hlcma0001/level2_marker_selection_table/",sample_type[s],"_SurfaceGenie_marker_gene_table_subset_arranged.rds"))
    
    table <- table %>% filter(table$annotation.l2 == cell_type)
    markers <- table$gene
    union_markers <- unique(sort(markers))
      
    
    union_markers <-unique(sort(union_markers))
    list <- list(union_markers)
    names(list)<- disease
    upset_list <- append(upset_list,list)
    
  }
  count <- table(unlist(upset_list))
  for_barplot <- data.frame(cbind(count))
  for_barplot$gene <- row.names(for_barplot)
  list2 <- list(for_barplot$gene)
  names(list2)<- cell_type
  wenn_list <- append(wenn_list,list2)
  
  #saveRDS(wenn_list,paste0(outdir, "/", "20220203_ourdataset_wenn_list_AT2.rds"))
  
  pdf(paste0(outdir,"/hlcma0001/level2_marker_selection_table/", cell_in_interesed[c], "_markers_upsetplot_barplot20220207.pdf"), width =8, height = 5)
  
  upset(fromList(upset_list),nintersects = NA,
        sets = sample_type,
        mainbar.y.label = "Sample Type Intersections", sets.x.label = "Markers per sample type", 
        text.scale = 1.2, mb.ratio = c(0.6, 0.4), order.by = "degree")
  
  ggplot(data=for_barplot, aes(x= reorder(gene,-count), y=count)) +
    geom_bar(stat="identity")+
    xlab(paste(length(for_barplot$gene),"marker genes"))+
    ylab("# of sample types")+
    ggtitle(paste0(cell_type," cell's markers across all sample types"))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 8),
          plot.title = element_text(hjust = 0.5))
  
  dev.off()
}


library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")[1:2]
library(RAM)
pdf(paste0(outdir, "/wenn_diagram_markersshart_20_exp.pdf"), width =10, height = 10)
group.venn(wenn_list, label=TRUE, 
           fill =myCol,
           cat.pos = c(-5, 10),
           lab.cex=0.8)
dev.off()

library(VennDiagram)
venn.diagram(
  x = wenn_list,
  category.names = names(wenn_list),
  filename = paste0(outdir, "/wenn_diagram_markersshart_20_exp.png"),
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600, 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  fontfamily = "serif",
  cex = 0.5,
  #categaory names
  cat.cex = 0.5,
  cat.default.pos = "text",
  cat.dist = c(0.04, -0.1),
  cat.fontfamily = "serif"
)



