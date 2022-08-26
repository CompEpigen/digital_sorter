
## Obtain all protein coding genes from chr 1-22 using biomaRt
## Reference: https://www.biostars.org/p/168203/

setwd("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter")
library(biomaRt)  
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
genes <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name","transcript_biotype"), 
                        filters = c("transcript_biotype","chromosome_name"),
                        values = list("protein_coding",c(1:22)), mart = mart)
prot.coding.genes <- unique(sort(genes$external_gene_name))

saveRDS(prot.coding.genes,"digital_sorter/data/prot.coding.genes.rds")
saveRDS(prot.coding.genes,"private_sorter/data/prot.coding.genes.rds")

