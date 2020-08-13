setwd("~/PycharmProjects/pepsickle-paper/data/raw")
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)

name_to_peptide_id_map <- getBM(attributes=c('external_gene_name', 'ensembl_peptide_id'), mart = ensembl)
name_to_peptide_id_map <- name_to_peptide_id_map[complete.cases(name_to_peptide_id_map),]
name_to_peptide_id_map <- name_to_peptide_id_map[name_to_peptide_id_map[,"ensembl_peptide_id"] != "",]
unique_id_map <- name_to_peptide_id_map[!duplicated(name_to_peptide_id_map[,'external_gene_name']),]

write.csv(unique_id_map, file="gene_protein_id_map.csv", row.names = F)
