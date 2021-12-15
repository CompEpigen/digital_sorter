library(dplyr)
library("viridis") 
library(ggplot2)
library(tidyr)
library(stringr)
library(pheatmap)
# Set working directory ####
wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)
rawdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/lung_mapped_cellxgene_fixed_addsplits"
datasetdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_each_dataset"
outdir = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/lung_mapped_cellxgene_fixed_addsplits/0_test/"


# Load functions
GetfileNames <- function(fileDir, pattern = ".rds"){
  filenames <- list.files(fileDir, pattern= pattern, full.names=TRUE)
  filename <- sub(paste0(".*",fileDir,"/"), "", filenames)
  filename <- sub(pattern, "", filename)  
  return(filename)
}

# Set parameters ####
list_samples_disease <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/LIST_lung_seuratsplit_7datasets_samples_disease.rds")
filenames <- GetfileNames(datasetdir)
cell_type <- c("Natural Killer","CD4+ Naive T","CD8+ Naive T","CD4+ Memory/Effector T",
               "CD8+ Memory/Effector T","B","Adventitial Fibroblast","Alveolar Fibroblast")
#d <- 3L
results <- data.frame()

for(d in 1:length(filenames)){
  print(filenames[d])
  sample_type <- list_samples_disease[[filenames[d]]]
  #s <- 1L #1-10
  for( s in 1:length(sample_type)){
    disease <- sample_type[s]
    print(disease)
    
    table <- readRDS(paste0(outdir,filenames[d],"/",sample_type[s],"_marker_gene_table_all.rds"))
    table <- table %>% filter(table$annotation.l2 %in% cell_type)
    table$group <- rep(paste0(filenames[d],":",sample_type[s]),length(row.names(table)))
    results <- rbind(results, table)
  }
  
  
}
saveRDS(results,paste0(outdir,"results_for_heatmap.RDS"))


select_col <- c("CD4+ Memory/Effector T CD3D","CD4+ Memory/Effector T CD4","CD4+ Naive T CD3D","CD4+ Naive T CD4",
                "CD8+ Naive T CD8B","CD8+ Naive T CD8A","CD8+ Naive T CD3D","B CD79A","Natural Killer NKG7",
                "CD8+ Memory/Effector T CD8B","CD8+ Memory/Effector T CD8A","CD8+ Memory/Effector T CD3D",
                "Adventitial Fibroblast ACTA2","Alveolar Fibroblast ACTA2")
col_order <- c("bischoff_2021:adjacent_normal","bischoff_2021:luad",
               "kim_2020:adjacent_normal" ,"kim_2020:nsclc_primary","kim_2020:nsclc_meta",
               "kim_2020:pleural_fluids","kim_2020:LymphNode_normal", "kim_2020:LymphNode_meta",
               "lambrechts_2018:adjacent_normal","lambrechts_2018:luad","lambrechts_2018:lusc" ,"lambrechts_2018:nsclc_primary",
               "travaglini_2020:general_normal","wu_2021:luad","wu_2021:lusc", "wu_2021:nsclc_primary","zilionis_2019:nsclc_primary" ,
               "zilionis_2019:blood", "song_2019:adjacent_normal" ,"song_2019:nsclc_primary")

## p-value ####
p_value <- results[,c(1,2,3,8)]
p_value$select <- paste(p_value$annotation.l2, p_value$gene)
p_value$select2 <- paste(p_value$group, p_value$gene)
p_value <- p_value %>% filter(p_value$select %in% select_col) %>% arrange(gene,group,desc(p_val))
p_value2 <- p_value[!duplicated(p_value$select2 ),]
p_value2 <- p_value2[,2:4]

p_value3 <- as.data.frame(pivot_wider(p_value2, names_from = group, values_from = p_val)) 
row.names(p_value3) <- p_value3$gene
p_value3$gene <- NULL
p_value4 <- p_value3[, col_order]
write.csv(p_value4,paste0(outdir,"p_value_results_for_heatmap.csv"))

## logFC ####
logFC <- results[,c(1,2,4,8)]
logFC$select <- paste(logFC$annotation.l2, logFC$gene)
logFC$select2 <- paste(logFC$group, logFC$gene)
logFC <- logFC %>% filter(logFC$select %in% select_col) %>% arrange(gene,group,desc(avg_log2FC))
logFC2 <- logFC[!duplicated(logFC$select2 ),]
logFC2 <- logFC2[,2:4]

logFC3 <- as.data.frame(pivot_wider(logFC2, names_from = group, values_from = avg_log2FC)) 
row.names(logFC3) <- logFC3$gene
logFC3$gene <- NULL
logFC4 <- logFC3[, col_order]
write.csv(logFC4,paste0(outdir,"logFC_results_for_heatmap.csv"))

## percentage ####
per <- results[,c(1,2,5,8)]
per$select <- paste(per$annotation.l2, per$gene)
per$select2 <- paste(per$group, per$gene)
per <- per %>% filter(per$select %in% select_col) %>% arrange(gene,group,desc(pct.1))
per2 <- per[!duplicated(per$select2 ),]
per2 <- per2[,2:4]

per3 <- as.data.frame(pivot_wider(per2, names_from = group, values_from = pct.1)) 
row.names(per3) <- per3$gene
per3$gene <- NULL
per4 <- per3[, col_order]
write.csv(per3,paste0(outdir,"percentage_results_for_heatmap.csv"))

df <- data.frame(x = colnames(per4))
annotation <- as.data.frame(str_split_fixed(df$x, ":", 2))
row.names(annotation) <- colnames(per4)
colnames(annotation) <- c("Cohort","Sample Type")
Var1 = rainbow(7)
names(Var1) = filenames
Var2 = c("lightblue", "darkred","navy","pink2","purple",
          "gold","orange","greenyellow","darkgreen","royalblue")
names(Var2) = unique(sort(annotation$"Sample Type"))

ann_colors = list(Cohort = Var1, "Sample Type" = Var2)

pdf(paste0(outdir,"Heatmap_cutoff.pdf"), width=10, height=6)

pheatmap::pheatmap(p_value4,
                   main = "p-value of well-known markers",
                   color = viridis(100),
                   display_numbers = TRUE, 
                   number_color = "white", 
                   number_format = '%.1e',
                   fontsize_number =5,
                   annotation = annotation,
                   annotation_colors = ann_colors, 
                   cluster_rows = F,
                   cluster_cols = F)
pheatmap::pheatmap(logFC4,
                   main = "log2FC of well-known markers",
                   color = viridis(100),
                   annotation_colors = ann_colors,
                   display_numbers = TRUE, 
                   number_color = "white", 
                   number_format = '%.2f',
                   annotation = annotation,
                   cluster_rows = F,
                   cluster_cols = F)
pheatmap::pheatmap(per4,
                   main = "percentage of well-known markers",
                   color = viridis(100),
                   display_numbers = TRUE, 
                   number_color = "white", 
                   number_format = '%.2f',
                   annotation = annotation,
                   annotation_colors = ann_colors,
                   cluster_rows = F,
                   cluster_cols = F)
dev.off()

pdf(paste0(outdir,"Heatmap_cutoff_cluster.pdf"), width=10, height=6)

pheatmap::pheatmap(p_value4,
                   main = "p-value of well-known markers",
                   color = viridis(100),
                   display_numbers = TRUE, 
                   number_color = "white", 
                   number_format = '%.1e', 
                   fontsize_number = 5,
                   annotation = annotation,
                   annotation_colors = ann_colors, 
                   cluster_rows = F,
                   cluster_cols = T)
pheatmap::pheatmap(logFC4,
                   main = "log2FC of well-known markers",
                   color = viridis(100),
                   annotation_colors = ann_colors,
                   display_numbers = TRUE, 
                   number_color = "white", 
                   number_format = '%.2f',
                   annotation = annotation,
                   cluster_rows = F,
                   cluster_cols = T)
pheatmap::pheatmap(per4,
                   main = "percentage of well-known markers",
                   color = viridis(100),
                   display_numbers = TRUE, 
                   number_color = "white", 
                   number_format = '%.2f',
                   annotation = annotation,
                   annotation_colors = ann_colors,
                   cluster_rows = F,
                   cluster_cols = T)
dev.off()
