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


