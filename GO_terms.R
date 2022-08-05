
library(dplyr)
library(biomaRt)

ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)


gene2GOtable <- read.table("GENE_LIST.txt", sep='', header = TRUE) 
 
gene2go <-getBM(attributes = c('hgnc_symbol','go_id','definition_1006','ensembl_gene_id', 'namespace_1003' ), values = gene2GOtable$ensembl_gene_id, filters = 'ensembl_gene_id', mart = ensembl)

gene2go_bio <- dplyr::filter(v,namespace_1003 == 'biological_process' )



