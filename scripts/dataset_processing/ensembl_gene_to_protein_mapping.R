#!/usr/bin/env RScript
"""
ensembl_gene_to_protein_mapping.R

For issues contact Ben Weeder (weeder@ohsu.edu)

This script pulls gene name and ensembl peptide id's for the human proteome. These are used downstream to map gene
names to unique searchable identifiers
"""
# TODO add relative file paths and I/O
# set working directory... this might need changed?
setwd("~/PycharmProjects/pepsickle-paper/data/raw")

library(biomaRt)
# use the main ensembl
ensembl <- useMart("ensembl")
# human genes
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)

# retrieve the common gene name and the unique ID
name_to_peptide_id_map <- getBM(attributes=c('external_gene_name', 'ensembl_peptide_id'), mart = ensembl)
# drop incomplete maps
name_to_peptide_id_map <- name_to_peptide_id_map[complete.cases(name_to_peptide_id_map),]
# drop mappings with no ID
name_to_peptide_id_map <- name_to_peptide_id_map[name_to_peptide_id_map[,"ensembl_peptide_id"] != "",]
#remove duplicates
unique_id_map <- name_to_peptide_id_map[!duplicated(name_to_peptide_id_map[,'external_gene_name']),]
# export
write.csv(unique_id_map, file="gene_protein_id_map.csv", row.names = F)
