library(Seurat)
library(SeuratDisk)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)

wd = "/Users/b370-mb03wl/Documents/P01_NSCLC/P01.2_cell_components/procdata"
setwd(wd)
ob = readRDS("lung_mapped_cellxgene.rds")
ob@meta.data$disease[!(ob@meta.data$disease=="adjacent normal"|ob@meta.data$disease=="normal"|ob@meta.data$disease=="nsclc")]="nsclc_meta"
table(ob@meta.data$disease)
table(ob@meta.data$donor)
ob@meta.data$dataset_origin = factor(ob@meta.data$dataset_origin, levels = c('song_2019','travaglini_2020','kim_2020'))
head(ob@meta.data$annotation.l1)
# CD45 (PTPRC)
features=c("PTPRC")
cols = c("steelblue","darkred","gold","coral2")
m1 =VlnPlot(ob, features, split.by = "disease", group.by = "annotation.l2", cols=cols,
        sort = TRUE,pt.size = 0, combine = FALSE)
m1
ggsave("CD45.pdf", height=6, width=15)
d = m1[[1]]$data
mu = aggregate(x=d$PTPRC,
          by=list(d$ident,d$split),
          FUN=mean, na.action=)
summary(mu)
mu = mu[mu$Group.2=='normal',]
set1 = unique(sort(mu[mu$x >= mean(mu$x, na.rm = T),]$Group.1))
set2 = levels(set1)[!levels(set1) %in% set1]
meta = ob$annotation.l2
meta[meta %in% set1] = paste0(features, "+")
meta[meta %in% set2] = paste0(features, "-")


ob <- AddMetaData(object = ob, metadata = meta, col.name = "split1")
ob_split <- SplitObject(ob, split.by = "disease")
length(ob_split)
for (i in 1:4){
  ob1 <- SplitObject(ob_split[[i]], split.by = "split1")
  features=c("CD8A","CD8B","CD4", "MRC1","CD14","SIGLEC1")
  a = DotPlot(ob1$`PTPRC+`, features = features, group.by = "annotation.l2") + RotatedAxis()
  b = VlnPlot(ob1$`PTPRC+`, features, stack=TRUE, flip=FALSE, cols = cols[i], 
          pt.size = 0, combine = FALSE, sort = FALSE, same.y.lims = TRUE,
          group.by = "annotation.l2", split.by = "disease") +
    theme(legend.position = "none")
  c = RidgePlot(ob1$`PTPRC+`, features, stack=TRUE, 
          combine = FALSE, sort = FALSE, same.y.lims = TRUE,
          group.by = "annotation.l2")
  pdf(paste0(names(ob_split)[i],"_","PTPRC+",".pdf"), width=15, height=12)
  grid.arrange(
    a, b, c,
    widths = c(1, 1),
    layout_matrix = rbind(c(1, 2),
                          c(3, 3)),
    top=textGrob(paste0(names(ob_split)[i], " (PTPRC+)"), gp=gpar(fontsize=20,font=3))
  )
  dev.off()

}

ob <- SplitObject(ob, split.by = "split1")
ob = ob$`PTPRC-`

features=c("EPCAM")
m1 =VlnPlot(ob, features, split.by = "disease", group.by = "annotation.l2", 
            cols = cols,sort = TRUE,pt.size = 0, combine = FALSE)
m1
ggsave("EPCAM.pdf", height=6, width=15)
d = m1[[1]]$data
mu = aggregate(x=d$EPCAM,
               by=list(d$ident,d$split),
               FUN=mean, na.action=)
summary(mu)
mu = mu[mu$Group.2=='normal',]
set1 = unique(sort(mu[mu$x >= mean(mu$x, na.rm = T),]$Group.1))
set2 = levels(set1)[!levels(set1) %in% set1]
meta = ob$annotation.l2
meta[meta %in% set1] = paste0(features, "+")
meta[meta %in% set2] = paste0(features, "-")

ob <- AddMetaData(object = ob, metadata = meta, col.name = "split2")
ob_split = SplitObject(ob, split.by = "disease")

features=c("NGFR","CD24","PDPN")
features=c("NGFR","CD24","PDPN","AGER","WIF1","ADGRF5","CEL","LTF","CEACAM6","TSPAN8","PLAT")
i=1
for (i in 1:4){
  print(names(ob_split)[i])
  ob1 <- SplitObject(ob_split[[i]], split.by = "split2")
  a = DotPlot(ob1$`EPCAM+`, features = features, group.by = "annotation.l2") + RotatedAxis()
  b = VlnPlot(ob1$`EPCAM+`, features, stack=TRUE, flip=FALSE, cols = cols[i], 
              pt.size = 0, combine = FALSE, sort = FALSE, same.y.lims=TRUE,
              group.by = "annotation.l2", split.by = "disease") +
    theme(legend.position = "none")
  c = RidgePlot(ob1$`EPCAM+`, features, stack=TRUE, 
                combine = FALSE, sort = FALSE, same.y.lims = TRUE,
                group.by = "annotation.l2")
    
  pdf(paste0(names(ob_split)[i],"_","EPCAM+",".pdf"), width=15, height=12)
  grid.arrange(
    a, b, c,
    widths = c(1, 1),
    layout_matrix = rbind(c(1, 2),
                          c(3, 3)),
    top=textGrob(paste0(names(ob_split)[i], " (EPCAM+)"), gp=gpar(fontsize=20,font=3))
  )
  dev.off()
  
}

i=4
for (i in 1:4){
  print(names(ob_split)[i])
  ob1 <- SplitObject(ob_split[[i]], split.by = "split2")
  seurat = ob1$`EPCAM+`
  # Feature selection
  seurat <- getMarkerGenes(
    seurat,
    assay = 'RNA',
    organism = 'hg',
    groups = c('annotation.l1'),
    name = 'EPCAM+_selection',
    only_pos = TRUE
  )
  names(seurat)
  summary(seurat)
  markers = seurat@misc$marker_genes$`EPCAM+_selection`$annotation.l1[seurat@misc$marker_genes$`EPCAM+_selection`$annotation.l1$on_cell_surface==TRUE,]
  write.csv(markers, paste0(names(ob_split)[i],"_","EPCAM+_selection.csv"))
  markers=markers[-1*grep("HLA",markers$gene),]
  features = rownames(markers[!duplicated(markers$annotation.l1),])
  names(features) = markers[!duplicated(markers$annotation.l1),]$annotation.l1
  features = features[order(names(features))]
  if (sum(features=="TFPI1")>0){ features[features=="TFPI1"] = "TFPI"}
  if (sum(features=="MIF1")>0){ features[features=="MIF1"] = "MIF"}
  if (sum(features=="LTF1")>0){ features[features=="LTF1"] = "LTF"}
  
  write.csv(features, paste0(names(ob_split)[i],"_","EPCAM+_selection_features4plot.csv"))
  features = unname(features)
  features = features[!duplicated(features)]
  
  a = DotPlot(ob1$`EPCAM+`, features = features, group.by = "annotation.l2") + RotatedAxis()
  b = VlnPlot(ob1$`EPCAM+`, features, stack=TRUE, flip=FALSE, cols = cols[i], 
              pt.size = 0, combine = FALSE, sort = FALSE, same.y.lims=TRUE,
              group.by = "annotation.l2", split.by = "disease") +
    theme(legend.position = "none")
  c = RidgePlot(ob1$`EPCAM+`, features, stack=TRUE, 
                combine = FALSE, sort = FALSE, same.y.lims = TRUE,
                group.by = "annotation.l2")
  
  pdf(paste0(names(ob_split)[i],"_","EPCAM+_selection",".pdf"), width=15, height=12)
  grid.arrange(
    a, b, c,
    widths = c(1, 1),
    layout_matrix = rbind(c(1, 2),
                          c(3, 3)),
    top=textGrob(paste0(names(ob_split)[i], " (EPCAM+)"), gp=gpar(fontsize=20,font=3))
  )
  dev.off()
  
}

## Level 3
ob <- SplitObject(ob, split.by = "split2")
ob = ob$`EPCAM-`

features=c("PECAM1")
m1 =VlnPlot(ob, features, split.by = "disease", group.by = "annotation.l2", 
            cols = c("steelblue","gold","darkred"),sort = TRUE,pt.size = 0, combine = FALSE)
m1
ggsave("PECAM1.pdf", height=6, width=15)
d = m1[[1]]$data
mu = aggregate(x=d$PECAM1,
               by=list(d$ident,d$split),
               FUN=mean, na.action=)
summary(mu)
mu = mu[mu$Group.2=='normal',]
set1 = unique(sort(mu[mu$x >= median(mu$x, na.rm = T),]$Group.1))
set2 = levels(set1)[!levels(set1) %in% set1]
meta = ob$annotation.l2
meta[meta %in% set1] = paste0(features, "+")
meta[meta %in% set2] = paste0(features, "-")

ob <- AddMetaData(object = ob, metadata = meta, col.name = "split3")
ob_split = SplitObject(ob, split.by = "disease")
cols = c("steelblue","darkred","gold")
features=c("CD34","CD44","PDPN","FLT4")

i=3
for (i in 1:3){
  print(names(ob_split)[i])
  ob1 <- SplitObject(ob_split[[i]], split.by = "split3")
  a = DotPlot(ob1$`PECAM1+`, features = features, group.by = "annotation.l2") + RotatedAxis()
  b = VlnPlot(ob1$`PECAM1+`, features, stack=TRUE, flip=FALSE, cols = cols[i], 
              pt.size = 0, combine = FALSE, sort = FALSE, same.y.lims=TRUE,
              group.by = "annotation.l2", split.by = "disease") +
    theme(legend.position = "none")
  c = RidgePlot(ob1$`PECAM1+`, features, stack=TRUE, 
                combine = FALSE, sort = FALSE, same.y.lims = TRUE,
                group.by = "annotation.l2")
  
  pdf(paste0(names(ob_split)[i],"_","PECAM1+",".pdf"), width=15, height=12)
  grid.arrange(
    a, b, c,
    widths = c(1, 1),
    layout_matrix = rbind(c(1, 2),
                          c(3, 3)),
    top=textGrob(paste0(names(ob_split)[i], " (PECAM1+)"), gp=gpar(fontsize=20,font=3))
  )
  dev.off()
  
}

features=c("CALD1","THY1","CD248","ITGA8","MCAM")
i=1
for (i in 1:3){
  print(names(ob_split)[i])
  ob1 <- SplitObject(ob_split[[i]], split.by = "split3")
  a = DotPlot(ob1$`PECAM1-`, features = features, group.by = "annotation.l2") + RotatedAxis()
  b = VlnPlot(ob1$`PECAM1-`, features, stack=TRUE, flip=FALSE, cols = cols[i], 
              pt.size = 0, combine = FALSE, sort = FALSE, same.y.lims=TRUE,
              group.by = "annotation.l2", split.by = "disease") +
    theme(legend.position = "none")
  c = RidgePlot(ob1$`PECAM1-`, features, stack=TRUE, 
                combine = FALSE, sort = FALSE, same.y.lims = TRUE,
                group.by = "annotation.l2")
  
  pdf(paste0(names(ob_split)[i],"_","PECAM1-",".pdf"), width=15, height=12)
  grid.arrange(
    a, b, c,
    widths = c(1, 1),
    layout_matrix = rbind(c(1, 2),
                          c(3, 3)),
    top=textGrob(paste0(names(ob_split)[i], " (PECAM1-)"), gp=gpar(fontsize=20,font=3))
  )
  dev.off()
  
}
