wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/dirty_scripts"
setwd(wd)

library(Seurat)
#library(SeuratDisk)
library(dplyr)
library(cerebroApp)
library(grid)
library(gridExtra)
library(ggplot2)

ob = readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/song_2019_addsplits.rds")
ob_split <- SplitObject(ob, split.by = "split1.1")

##1
ob_split_PTPRC = ob_split$`PTPRC+`
seurat <- getMarkerGenes(
  ob_split_PTPRC,assay = 'RNA',
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
  marker_gene_table <- marker_gene_table %>% filter(marker_gene_table$on_cell_surface =="TRUE")%>%
    arrange(desc(avg_log2FC))
  
saveRDS(marker_gene_table,"marker_level1_CD45pos.rds")

table(marker_gene_table$annotation.l2)
length(table(marker_gene_table$annotation.l2)) #18
#16 cell types has marker
repeat_genes <- as.data.frame(table(marker_gene_table$gene) <3 )
colnames(repeat_genes) <- "repeat3"
repeat_genes$gene <- row.names(repeat_genes)
marker_gene_table <- merge(marker_gene_table,repeat_genes,by = "gene", all.x = T) 
marker_gene_table <- marker_gene_table[marker_gene_table$repeat3,]
marker_gene_table <- marker_gene_table %>% 
  arrange(annotation.l2,desc(avg_log2FC)) 
marker_gene_table <- marker_gene_table[!duplicated(marker_gene_table$"annotation.l2"),]%>% 
  arrange(desc(avg_log2FC)) 
saveRDS(marker_gene_table,"marker_level1_CD45pos_selected.rds")

##2
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
  test = 'wilcox',
  verbose = TRUE
)
#marker_gene <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]][["gene"]]
marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
marker_gene_table <- marker_gene_table %>% filter(marker_gene_table$on_cell_surface =="TRUE")%>%
  arrange(desc(avg_log2FC))

saveRDS(marker_gene_table,"marker_level2_CD45neg.rds")
table(marker_gene_table$annotation.l2)
length(table(marker_gene_table$annotation.l2)) #8
#16 cell types has marker
table(marker_gene_table$gene)
repeat_genes <- as.data.frame(table(marker_gene_table$gene) <3 )
colnames(repeat_genes) <- "repeat3"
repeat_genes$gene <- row.names(repeat_genes)
marker_gene_table <- merge(marker_gene_table,repeat_genes,by = "gene", all.x = T) 
marker_gene_table <- marker_gene_table[marker_gene_table$repeat3,]
marker_gene_table <- marker_gene_table %>% 
  arrange(annotation.l2,desc(avg_log2FC)) 
marker_gene_table <- marker_gene_table[!duplicated(marker_gene_table$"annotation.l2"),]%>% 
  arrange(desc(avg_log2FC)) 
saveRDS(marker_gene_table,"marker_level2_CD45neg_selected.rds")

##3
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
marker_gene_table <- marker_gene_table %>% filter(marker_gene_table$on_cell_surface =="TRUE")%>%
  arrange(desc(avg_log2FC))

saveRDS(marker_gene_table,"marker_level3_EPCAMneg.rds")
table(marker_gene_table$annotation.l2)
length(table(marker_gene_table$annotation.l2)) #8
#16 cell types has marker
table(marker_gene_table$gene)
repeat_genes <- as.data.frame(table(marker_gene_table$gene) <3 )
colnames(repeat_genes) <- "repeat3"
repeat_genes$gene <- row.names(repeat_genes)
marker_gene_table <- merge(marker_gene_table,repeat_genes,by = "gene", all.x = T) 
marker_gene_table <- marker_gene_table[marker_gene_table$repeat3,]
marker_gene_table <- marker_gene_table %>% 
  arrange(annotation.l2,desc(avg_log2FC)) 
marker_gene_table <- marker_gene_table[!duplicated(marker_gene_table$"annotation.l2"),]%>% 
  arrange(desc(avg_log2FC)) 
saveRDS(marker_gene_table,"marker_level3_EPCAMneg_selected.rds")
##4
ob4 <- SplitObject(ob, split.by = "split3")
ob_split_PECAM1 = ob4$`PECAM1-`
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
marker_gene_table <- marker_gene_table %>% filter(marker_gene_table$on_cell_surface =="TRUE")%>%
  arrange(desc(avg_log2FC))

saveRDS(marker_gene_table,"marker_level4_PECAM1neg.rds")
table(marker_gene_table$annotation.l2)
length(table(marker_gene_table$annotation.l2)) #20
#18 cell types has marker
table(marker_gene_table$gene)
repeat_genes <- as.data.frame(table(marker_gene_table$gene) <3 )
colnames(repeat_genes) <- "repeat3"
repeat_genes$gene <- row.names(repeat_genes)
marker_gene_table <- merge(marker_gene_table,repeat_genes,by = "gene", all.x = T) 
marker_gene_table <- marker_gene_table[marker_gene_table$repeat3,]
marker_gene_table <- marker_gene_table %>% 
  arrange(annotation.l2,desc(avg_log2FC)) 
marker_gene_table <- marker_gene_table[!duplicated(marker_gene_table$"annotation.l2"),]%>% 
  arrange(desc(avg_log2FC)) 
saveRDS(marker_gene_table,"marker_level4_PECAM1neg_selected.rds")

if(F){ 
## extract top 10 gene depends on logFC ####
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