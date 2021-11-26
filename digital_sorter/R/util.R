##--------------------------------------------------------------------------##
## Get list of genes in cell surface through gene ontology term GO:0009986.
##--------------------------------------------------------------------------##

get_cell_markers = function(dataset="hsapiens_gene_ensembl", GO = 'GO:0009986'){
  
  ensembl = biomaRt::useEnsembl("ensembl",dataset= "hsapiens_gene_ensembl") #uses human ensembl annotations
  #gets gene symbol, transcript_id and go_id for all genes annotated with GO:0009986
  gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'go_id'),
                     filters = 'go_parent_term', values = GO , mart = ensembl)
  
  gene.data.unique <- gene.data[!duplicated(gene.data$hgnc_symbol),]
  cell.surface.marker = gene.data.unique$hgnc_symbol
}

if(F){ 
  #grap cell ontology file
  library(stringr)
  library(dplyr)
  setwd("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter")
  CL<- read.csv("CL.csv")
  CL<- CL[,1:2]
  ID <- as.character(CL$Class.ID)
  ID2 <- as.data.frame(strsplit(ID,'_'))
  ID3 <- data.table::transpose(ID2)[,1]
  ID4 <- str_trunc(ID3,2,side="left",ellipsis="")
  CL$class <- ID4
  CL_CL <- CL %>% filter(CL$class == "CL")
  write.csv(CL_CL,"procdata/CL_cell_names_lee.csv")
}


if(F){ 
  ## Read RDS files ####
  ReadRDSFiles <- function(fileDir, envir = .GlobalEnv) {
    
    # pattern for file searching
    p <- ".rds$"
    rds <- list.files(path = fileDir, pattern = p)
    out <- vapply(rds, FUN = function(.x) {
      nm <- sub(pattern = p, replacement = "", x = .x)
      # read data in
      # out <- readRDS(file = x)
      
      # load in global env
      assign(nm, value = readRDS(file = file.path(fileDir, .x)), envir = envir)
      if (!exists(nm, envir = envir)) return(FALSE)
      TRUE
    }, FUN.VALUE = logical(1L), USE.NAMES = FALSE)
    
    if (!all(out)) warning("Some `.rds` files not loaded.", call. = FALSE)
    
    spc <- paste0(rep('*', times = nchar(fileDir) + 1), collapse = "")
    cat("RDS Files loaded from:", fileDir, spc, rds[out], sep = "\n ", spc)
  }
  ## Get files names####
  GetfileNames <- function(fileDir, pattern = ".rds"){
    filenames <- list.files(fileDir, pattern= pattern, full.names=TRUE)
    filename <- sub(paste0(".*",fileDir,"/"), "", filenames)
    filename <- sub(pattern, "", filename)  
    return(filename)
  }
}