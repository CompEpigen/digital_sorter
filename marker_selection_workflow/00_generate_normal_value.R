library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)

wd = "/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter"
setwd(wd)
ob = readRDS("/omics/groups/OE0219/internal/MJMC/P01_NSCLC/P01.2_enriched_cell_components/analysis/azimuth-meta-analysis/human_lung/cellxgene/lung_mapped_cellxgene_fixed.rds")
table(ob@meta.data$disease)
table(ob@meta.data$donor)
table(ob@meta.data$dataset_origin)

# CD45 (PTPRC)
features=c("PTPRC")
cols = c("steelblue","darkred","gold","coral2","pink","yellowgreen","violet","blue","grey","black")
m1 =VlnPlot(ob, features, split.by = "disease", group.by = "annotation.l2", cols=cols,
        sort = TRUE,pt.size = 0, combine = FALSE)
m1
d = m1[[1]]$data
mu = aggregate(x=d$PTPRC,
          by=list(d$ident,d$split),
          FUN=mean)
summary(mu)
mu = mu[mu$Group.2=='general_normal',] #categorize PTPRC+/-
saveRDS(mu,"dirty_scripts/normal_mean_PTPRC.rds")
set1 = unique(sort(mu[mu$x >= mean(mu$x, na.rm = T),]$Group.1))
set2 = levels(set1)[!levels(set1) %in% set1]
meta = ob$annotation.l2
meta[meta %in% set1] = paste0(features, "+")
meta[meta %in% set2] = paste0(features, "-")


ob <- AddMetaData(object = ob, metadata = meta, col.name = "split1")

ob <- SplitObject(ob, split.by = "split1")
ob = ob$`PTPRC-`

features=c("EPCAM")
m1 =VlnPlot(ob, features, split.by = "disease", group.by = "annotation.l2", 
            cols = cols,sort = TRUE,pt.size = 0, combine = FALSE)
m1

d = m1[[1]]$data
mu = aggregate(x=d$EPCAM,
               by=list(d$ident,d$split),
               FUN=mean)
summary(mu)
mu = mu[mu$Group.2=='general_normal',]
saveRDS(mu,"dirty_scripts/normal_mean_EPCAM.rds")
set1 = unique(sort(mu[mu$x >= mean(mu$x, na.rm = T),]$Group.1))
set2 = levels(set1)[!levels(set1) %in% set1]
meta = ob$annotation.l2
meta[meta %in% set1] = paste0(features, "+")
meta[meta %in% set2] = paste0(features, "-")

ob <- AddMetaData(object = ob, metadata = meta, col.name = "split2")

## Level 3
ob <- SplitObject(ob, split.by = "split2")
ob = ob$`EPCAM-`

features=c("PECAM1")
m1 =VlnPlot(ob, features, split.by = "disease", group.by = "annotation.l2", 
            cols = cols,sort = TRUE,pt.size = 0, combine = FALSE)
m1

d = m1[[1]]$data
mu = aggregate(x=d$PECAM1,
               by=list(d$ident,d$split),
               FUN=mean)
summary(mu)
mu = mu[mu$Group.2=='general_normal',]
saveRDS(mu,"dirty_scripts/normal_mean_PECAM1.rds")
set1 = unique(sort(mu[mu$x >= mean(mu$x, na.rm = T),]$Group.1))
set2 = levels(set1)[!levels(set1) %in% set1]
meta = ob$annotation.l2
meta[meta %in% set1] = paste0(features, "+")
meta[meta %in% set2] = paste0(features, "-")

ob <- AddMetaData(object = ob, metadata = meta, col.name = "split3")

## Level 4
ob <- SplitObject(ob, split.by = "split3")
ob = ob$`PECAM1-`

features=c("MME")
m1 =VlnPlot(ob, features, split.by = "disease", group.by = "annotation.l2", 
            cols = cols ,sort = TRUE,pt.size = 0, combine = FALSE)
m1

d = m1[[1]]$data
mu = aggregate(x=d$MME,
               by=list(d$ident,d$split),
               FUN=mean)
summary(mu)
mu = mu[mu$Group.2=='general_normal',]
saveRDS(mu,"dirty_scripts/normal_mean_MME.rds")
set1 = unique(sort(mu[mu$x >= median(mu$x, na.rm = T),]$Group.1))
set2 = levels(set1)[!levels(set1) %in% set1]
meta = ob$annotation.l2
meta[meta %in% set1] = paste0(features, "+")
meta[meta %in% set2] = paste0(features, "-")
ob <- AddMetaData(object = ob, metadata = meta, col.name = "split4")
