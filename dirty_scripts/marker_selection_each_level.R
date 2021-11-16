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
    only_pos = TRUE,
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
  
saveRDS(marker_gene_table,"marker_CD45pos_level1.rds")


##2
ob2 <- SplitObject(ob, split.by = "split2")
ob_split_PTPRC0 = ob2$`EPCAM+`

seurat <- getMarkerGenes(
  ob_split_PTPRC0,assay = 'RNA',
  organism = 'hg',
  groups = c('annotation.l2'),
  name = 'cerebro_seurat',
  only_pos = TRUE,
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

saveRDS(marker_gene_table,"marker_CD45neg_level2.rds")

##3
ob3 <- SplitObject(ob, split.by = "split3")
ob_split_EPCAM0 = ob3$`PECAM1+` #EPCAM-
seurat <- getMarkerGenes(
  ob_split_EPCAM0,assay = 'RNA',
  organism = 'hg',
  groups = c('annotation.l2'),
  name = 'cerebro_seurat',
  only_pos = TRUE,
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

saveRDS(marker_gene_table,"marker_EPCAMneg_level3.rds")

##4
ob4 <- SplitObject(ob, split.by = "split3")
ob_split_PECAM1 = ob4$`PECAM1-`
seurat <- getMarkerGenes(
  ob_split_PECAM1,assay = 'RNA',
  organism = 'hg',
  groups = c('annotation.l2'),
  name = 'cerebro_seurat',
  only_pos = TRUE,
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

saveRDS(marker_gene_table,"marker_PECAM1neg_level4.rds")

if(F){ 
## extract top 10 gene depends on logFC ####
level1 <- readRDS("marker_CD45pos_level1.rds")
gene_dge <- unique(level1$gene)[1:10]
#by logFC ("AREG","CLU","HLA-DPB1","HLA-DRB1","HLA-DPA1","LGALS3","HLA-DRA","CD74", "PLAUR","MRC1") 
###by subset:("AREG","HLA-DPB1","HLA-DRB1","HLA-DPA1","HLA-DRA","CD74","PLAUR","MRC1", "TYROBP","FCER1G-E") 
level1 <- level1 %>% arrange(p_val)
gene_dge <- unique(level1$gene)[1:10]
#by p-value ("LGALS3","MRC1","CTSZ","CD9","TYROBP","HLA-DMA","LAMP1,"CD63", "HLA-DRA","ANXA2") 
DotPlot(ob_split_PTPRC, features = gene_dge, group.by = "annotation.l2") + RotatedAxis()+ theme(legend.text=element_text(size=12),
                                                                                       axis.text=element_text(size=12),
                                                                                       axis.title=element_text(size=14),
                                                                                       legend.title=element_text(size=12))
level2 <- readRDS("marker_CD45neg_level2.rds")
gene_dge <- unique(level2$gene)[1:10]
#by logFC ("SPARC","HLA-E","LY6D","LGALS1","HLA-B","B2M","EGFL7","VAMP5", "HYAL2","CLU") 
###by subset:c("SPARC","EGFL7","HYAL2","VAMP5","CLU","C3","HLA-E", "BST2","CD24","TFPI") 
level2 <- level2 %>% arrange(p_val)
gene_dge <- unique(level2$gene)[1:10]
#by p-value ("TSPAN8","CD59","ADGRF5","GPC3","CD24","CEACAM6","LGALS3,"SPARC", "ANXA1","B2M") 
DotPlot(ob_split_PTPRC0, features = gene_dge, group.by = "annotation.l2") + RotatedAxis()+ theme(legend.text=element_text(size=12),
                                                                                                axis.text=element_text(size=12),
                                                                                                axis.title=element_text(size=14),
                                                                                                legend.title=element_text(size=12))


level3 <- readRDS("marker_EPCAMneg_level3.rds")
gene_dge <- unique(level3$gene)[1:10]
#by logFC ("CLU","TFPI","SPARC","CD9","TGFBR2","SDC2","HLA-DRB1","ENG", "ICAM1","ANXA2") 
###by subset:c("SPARC","CLU","TFPI","CD59","HYAL2","TGFBR2", "APP", "ENG","LGALS1" ,"CD9") 
level3 <- level3 %>% arrange(p_val)
gene_dge <- unique(level3$gene)[1:10]
#by p-value ("MRC1","ENG","TFPI","SPARC","CLU","TGFBR2","SDC2,"CD9", "ICAM1","HLA-DRB1") 
DotPlot(ob_split_EPCAM0, features = gene_dge, group.by = "annotation.l2") + RotatedAxis()+ theme(legend.text=element_text(size=12),
                                                                                                 axis.text=element_text(size=12),
                                                                                                 axis.title=element_text(size=14),
                                                                                                 legend.title=element_text(size=12))

level4 <- readRDS("marker_PECAM1neg_level4.rds")
gene_dge <- unique(level4$gene)[1:10]
#by logFC ("SPARC","ITGB1","SDC2","LGALS1","CD59","CD74","TFPI","CD63", "HSP90AB1","HMGB1") 
###by subset: c("SPARC","SDC2","CD59","LGALS1","TFPI","ITGB1","CD63","CALR", "BST2","TIMP2") 
level4 <- level4 %>% arrange(p_val)
gene_dge <- unique(level4$gene)[1:10]
#by p-value ("MRC1","ENG","TFPI","SPARC","CLU","TGFBR2","SDC2,"CD9", "ICAM1","HLA-DRB1") 
DotPlot(ob_split_PECAM1, features = gene_dge, group.by = "annotation.l2") + RotatedAxis()+ theme(legend.text=element_text(size=12),
                                                                                                 axis.text=element_text(size=12),
                                                                                                 axis.title=element_text(size=14),
                                                                                                 legend.title=element_text(size=12))
}