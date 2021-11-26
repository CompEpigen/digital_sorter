wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/dirty_scripts"
setwd(wd)

library(Seurat)
library(dplyr)
library(cerebroApp)
library(grid)
library(gridExtra)
library(ggplot2)
library(DESeq2)

ob = readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/song_2019_addsplits.rds")
ob_split <- SplitObject(ob, split.by = "split1.1")

## Level 1 ####
ob_split_PTPRC = ob_split$`PTPRC+`

seurat <- cerebroApp::getMarkerGenes(
  ob_split_PTPRC,assay = 'RNA',
    organism = 'hg',
    groups = c('annotation.l2'),
    name = 'cerebro_seurat',
    only_pos = F,
    min_pct = 0.7,
    thresh_logFC = 0.25,
    thresh_p_val = 0.01,
    test = 'wilcox', #DESeq2 no markers found
    verbose = TRUE
  )
marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]] #3708
marker_gene_table <- marker_gene_table %>% 
  filter(marker_gene_table$on_cell_surface =="TRUE")#280

# step 0: filter for FDR < 0.05
marker_gene_table <- marker_gene_table %>%
  filter(marker_gene_table$p_val_adj<0.05)%>% #224
  arrange(desc(avg_log2FC))
  
saveRDS(marker_gene_table,"marker_level1_CD45pos.rds")

features <- unique(sort(marker_gene_table$gene))

# plot the initial dot plot
pdf("Dotplot_level1.pdf", width=14, height=6)
d <- DotPlot(ob_split_PTPRC, features = features, group.by = "annotation.l2") + RotatedAxis()
d
dev.off()

# step 1: filter for the average expression > 1
exp_data <- d[["data"]]
exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
exp_data <- exp_data[,c(1,2,5,6)]
marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")

marker_gene_table2 <- marker_gene_table2 %>%
  filter(marker_gene_table2$avg.exp >1)   #222

# step 1.5: filter out frequently appearing markers
table(marker_gene_table2$annotation.l2)
length(table(marker_gene_table2$annotation.l2)) #18
table(marker_gene_table2$gene)
repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <=5 )
colnames(repeat_genes) <- "repeat5"
repeat_genes$gene <- row.names(repeat_genes)
marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat5,] #51

# step 2: arrange according to abs(log2FC)
marker_gene_table3 <- marker_gene_table3 %>% 
  arrange(annotation.l2,desc(abs(avg_log2FC)))

marker_gene_table3$mergeID <- NULL


# step 3: select top 10 for each cell type (exclude HLA related)
library(data.table)
subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #51
subset_marker_gene=subset_marker_gene[-1*grep("HLA",subset_marker_gene$gene),] #41
subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2]

saveRDS(subset_marker_gene,"marker_level1_CD45pos_selected.rds")

## Level 2 ####
ob = ob_split$`PTPRC-`
ob2 <- SplitObject(ob, split.by = "split2")
ob_split_PTPRC0 = ob2$`EPCAM+`

seurat <- getMarkerGenes(
  ob_split_PTPRC0,assay = 'RNA',
  organism = 'hg',
  groups = c('annotation.l2'),
  name = 'cerebro_seurat',
  only_pos = F,
  min_pct = 0.7,
  thresh_logFC = 0.25,
  thresh_p_val = 0.01,
  test = 'wilcox', #DESeq2 no markers found
  verbose = TRUE
)
marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
marker_gene_table <- marker_gene_table %>% 
  filter(marker_gene_table$on_cell_surface =="TRUE") #267

# step 0: filter for FDR < 0.05
marker_gene_table <- marker_gene_table %>%
  filter(marker_gene_table$p_val_adj<0.05)%>% #201
  arrange(desc(avg_log2FC))

saveRDS(marker_gene_table,"marker_level2_CD45neg.rds")
features <- unique(sort(marker_gene_table$gene))

# plot the initial dot plot
pdf("Dotplot_level2.pdf", width=20, height=6)
d <- DotPlot(ob_split_PTPRC0, features = features, group.by = "annotation.l2") + RotatedAxis()
d
dev.off()

# step 1: filter for the average expression > 1
exp_data <- d[["data"]]
exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
exp_data <- exp_data[,c(1,2,5,6)]
marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")

marker_gene_table2 <- marker_gene_table2 %>%
  filter(marker_gene_table2$avg.exp >1)   #182

# step 1.5: filter out frequently appearing markers
table(marker_gene_table2$annotation.l2)
length(table(marker_gene_table2$annotation.l2)) #11
table(marker_gene_table2$gene)
repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <=5 )
colnames(repeat_genes) <- "repeat5"
repeat_genes$gene <- row.names(repeat_genes)
marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat5,] #150

# step 2: arrange according to abs(log2FC)
marker_gene_table3 <- marker_gene_table3 %>% 
  arrange(annotation.l2,desc(abs(avg_log2FC)))

marker_gene_table3$mergeID <- NULL


# step 3: select top 10 for each cell type
library(data.table)
subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #150
subset_marker_gene=subset_marker_gene[-1*grep("HLA",subset_marker_gene$gene),] #127
subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] #86

saveRDS(subset_marker_gene,"marker_level2_CD45neg_selected.rds")


## Level 3 ####
ob = ob2$`EPCAM-`
ob3 <- SplitObject(ob, split.by = "split3")
ob_split_EPCAM0 = ob3$`PECAM1+` #EPCAM-
seurat <- getMarkerGenes(
  ob_split_EPCAM0,assay = 'RNA',
  organism = 'hg',
  groups = c('annotation.l2'),
  name = 'cerebro_seurat',
  only_pos = F,
  min_pct = 0.7,
  thresh_logFC = 0.25,
  thresh_p_val = 0.01,
  test = 'wilcox', 
  verbose = TRUE
)
marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
marker_gene_table <- marker_gene_table %>% 
  filter(marker_gene_table$on_cell_surface =="TRUE") #71

# step 0: filter for FDR < 0.05
marker_gene_table <- marker_gene_table %>%
  #filter(marker_gene_table$p_val_adj<0.1)%>% 
  arrange(desc(avg_log2FC))

saveRDS(marker_gene_table,"marker_level3_EPCAMneg_noFDR.rds")

features <- unique(sort(marker_gene_table$gene))
# plot the initial dot plot
pdf("Dotplot_level3_noFDR.pdf", width=11, height=6)
d <- DotPlot(ob_split_EPCAM0, features = features, group.by = "annotation.l2") + RotatedAxis()
d
dev.off()

# step 1: filter for the average expression > 1
exp_data <- d[["data"]]
exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
exp_data <- exp_data[,c(1,2,5,6)]
marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")

marker_gene_table2 <- marker_gene_table2 %>%
  filter(marker_gene_table2$avg.exp >1)   #65

# step 1.5: filter out frequently appearing markers
table(marker_gene_table2$annotation.l2)
length(table(marker_gene_table2$annotation.l2)) #7
table(marker_gene_table2$gene)
repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <=5 )
colnames(repeat_genes) <- "repeat5"
repeat_genes$gene <- row.names(repeat_genes)
marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat5,] #65

# step 2: arrange according to abs(log2FC)
marker_gene_table3 <- marker_gene_table3 %>% 
  arrange(annotation.l2,desc(abs(avg_log2FC)))

marker_gene_table3$mergeID <- NULL


# step 3: select top 10 for each cell type
library(data.table)
subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #65
subset_marker_gene=subset_marker_gene[-1*grep("HLA",subset_marker_gene$gene),] #48
subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] #45

saveRDS(marker_gene_table,"marker_level3_EPCAMneg_selected_noFDR.rds")

## Level 4 ####
ob_split_PECAM1 = ob3$`PECAM1-`
ob4 <- SplitObject(ob_split_PECAM1, split.by = "split4")
seurat <- getMarkerGenes(
  ob_split_PECAM1,assay = 'RNA',
  organism = 'hg',
  groups = c('annotation.l2'),
  name = 'cerebro_seurat',
  only_pos = F,
  min_pct = 0.7,
  thresh_logFC = 0.25,
  thresh_p_val = 0.01,
  test = 'wilcox',
  verbose = TRUE
)
marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
marker_gene_table <- marker_gene_table %>% 
  filter(marker_gene_table$on_cell_surface =="TRUE") #76

# step 0: filter for FDR < 0.05
marker_gene_table <- marker_gene_table %>%
  filter(marker_gene_table$p_val_adj<0.05)%>% #41
  arrange(desc(avg_log2FC))

saveRDS(marker_gene_table,"marker_level4_PECAM1neg.rds")
features <- unique(sort(marker_gene_table$gene))

# plot the initial dot plot
pdf("Dotplot_level4.pdf", width=18, height=6)
d <- DotPlot(ob_split_PECAM1, features = features, group.by = "annotation.l2") + RotatedAxis()
d
dev.off()

# step 1: filter for the average expression > 1
exp_data <- d[["data"]]
exp_data$mergeID <- paste(exp_data$features.plot,exp_data$id)
exp_data <- exp_data[,c(1,2,5,6)]
marker_gene_table$mergeID <- paste(marker_gene_table$gene, marker_gene_table$annotation.l2)
marker_gene_table2 <- merge(marker_gene_table, exp_data, by = "mergeID")

marker_gene_table2 <- marker_gene_table2 %>%
  filter(marker_gene_table2$avg.exp >1)   #29

# step 1.5: filter out frequently appearing markers
table(marker_gene_table2$annotation.l2)
length(table(marker_gene_table2$annotation.l2)) #3
table(marker_gene_table2$gene)
repeat_genes <- as.data.frame(table(marker_gene_table2$gene) <=5 )
colnames(repeat_genes) <- "repeat5"
repeat_genes$gene <- row.names(repeat_genes)
marker_gene_table3 <- merge(marker_gene_table2,repeat_genes,by = "gene", all.x = T) 
marker_gene_table3 <- marker_gene_table3[marker_gene_table3$repeat5,] #29

# step 2: arrange according to abs(log2FC)
marker_gene_table3 <- marker_gene_table3 %>% 
  arrange(annotation.l2,desc(abs(avg_log2FC)))

marker_gene_table3$mergeID <- NULL


# step 3: select top 10 for each cell type
library(data.table)
subset_marker_gene <- data.table(marker_gene_table3, key="annotation.l2") #29
subset_marker_gene=subset_marker_gene[-1*grep("HLA",subset_marker_gene$gene),] #28
subset_marker_gene <- subset_marker_gene[, head(.SD, 10), by=annotation.l2] #20
saveRDS(marker_gene_table,"marker_level4_PECAM1neg_selected.rds")


## extract top 10 gene depends on logFC ####
if(F){ 

level1 <- readRDS("marker_level1_CD45pos_selected.rds")
gene_dge <- unique(level1$gene)[1:10]
#by logFC
###by subset
level1 <- level1 %>% arrange(p_val)
gene_dge <- unique(level1$gene)[1:10]
#by p-value 
DotPlot(ob_split_PTPRC, features = gene_dge, group.by = "annotation.l2") + RotatedAxis()+ theme(legend.text=element_text(size=12),
                                                                                       axis.text=element_text(size=12),
                                                                                       axis.title=element_text(size=14),
                                                                                       legend.title=element_text(size=12))
level2 <- readRDS("marker_level2_CD45neg_selected.rds")
gene_dge <- unique(level2$gene)[1:10]
#by logFC 
##by subset:
level2 <- level2 %>% arrange(p_val)
gene_dge <- unique(level2$gene)[1:10]
#by p-value
DotPlot(ob_split_PTPRC0, features = gene_dge, group.by = "annotation.l2") + RotatedAxis()+ theme(legend.text=element_text(size=12),
                                                                                                axis.text=element_text(size=12),
                                                                                                axis.title=element_text(size=14),
                                                                                                legend.title=element_text(size=12))


level3 <- readRDS("marker_level3_EPCAMneg_selected.rds")
gene_dge <- unique(level3$gene)[1:10]
#by logFC
###by subset
level3 <- level3 %>% arrange(p_val)
gene_dge <- unique(level3$gene)[1:10]
#by p-value
DotPlot(ob_split_EPCAM0, features = gene_dge, group.by = "annotation.l2") + RotatedAxis()+ theme(legend.text=element_text(size=12),
                                                                                                 axis.text=element_text(size=12),
                                                                                                 axis.title=element_text(size=14),
                                                                                                 legend.title=element_text(size=12))

level4 <- readRDS("marker_level4_PECAM1neg_selected.rds")
gene_dge <- unique(level4$gene)[1:10]
#by logFC
###by subset
level4 <- level4 %>% arrange(p_val)
gene_dge <- unique(level4$gene)[1:10]
#by p-value 
DotPlot(ob_split_PECAM1, features = gene_dge, group.by = "annotation.l2") + RotatedAxis()+ theme(legend.text=element_text(size=12),
                                                                                                 axis.text=element_text(size=12),
                                                                                                 axis.title=element_text(size=14),
                                                                                                 legend.title=element_text(size=12))
}