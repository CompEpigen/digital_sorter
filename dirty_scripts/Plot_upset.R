

library(Seurat)
library(dplyr)
library(data.table)

wd = "/omics/groups/OE0219/internal/MJMC/P01_NSCLC/P01.2_enriched_cell_components/analysis/azimuth-meta-analysis/human_lung/results/by_cohorts_addsplits"
setwd(wd)

split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))


system2(command = "ls", args = paste0("*/*/marker_selected_all_levels.rds > ", "tmp.file"))
in.files = read.table("tmp.file")
in.files
i=4
ctypes = c()
bdf = data.frame()
for (i in 1:nrow(in.files)){
  f = in.files[i,]  
  df = as.data.frame(readRDS(f))
  class(df)
  ctypes = c(ctypes, as.vector(df$annotation.l2))
  folders = split_path(f)[2:3]
  bdf = rbind(bdf, data.frame(df, folders[1], folders[2]))
}
ctypes
ctypes = unique(sort(ctypes))
head(bdf)
tail(bdf)

library(UpSetR)
ctypes
i = 4
for (i in 1:length(ctypes)){
  ctype = ctypes[i]  
  print(ctype)
  df = bdf[bdf$annotation.l2==ctype,]
  tdf = df[,c("folders.1.","gene")]
  xy.df = tdf[,2]
  xy.list
  xy.list <- split(xy.df, factor(tdf$folders.1.))
  
  #-----place holder for a function
  #need to sort and unique the gene list
  #-----
  
  upset(fromList(xy.list), order.by = "freq", nintersects = 1000)
  table(unlist(xy.list))
  
  #-----place holder for a function
  # barplot: x-axis is genes; y-axis number of sets/groups
  #-----
}
