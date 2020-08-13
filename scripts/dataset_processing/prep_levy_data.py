#!usr/bin/env python3
"""

"""

import pandas as pd
import re
from extraction_functions import *
from optparse import OptionParser


# define command line parameters
parser = OptionParser()
parser.add_option("-i", "--in-file",
                  help="CSV of data from levy et al.")
parser.add_option("-o", "--out",
                  help="output directory where processed CSV is exported")
parser.add_option("-a", "--annotation_file", dest="annotation_file",
                  help="CSV of ")

(options, args) = parser.parse_args()

# load in data

dat = pd.read_csv(options.in_file, low_memory=False, header=2, skiprows=[3],
                  usecols=range(0, 4))

annot = pd.read_csv(options.annotation_file, low_memory=False)

# fix column names
dat.columns = ['fragment', 'gene_names', 'logF_med_intensity_UT',
               'logF_med_intensity_T']

# remove decoy peptides
dat = dat[(dat['logF_med_intensity_T'].notna()) |
          (dat['logF_med_intensity_UT'].notna())]

levy_df = pd.DataFrame(columns=['fragment', 'gene_names', 'intensity',
                                'Proteasome'])

for row in dat.iterrows():
    if row[1]['logF_med_intensity_UT'] > 0:
        new_row = pd.Series([row[1]['fragment'], row[1]['gene_names'],
                             row[1]['logF_med_intensity_UT'], 'C'],
                            index=levy_df.columns)
        levy_df = levy_df.append(new_row, ignore_index=True)
    if row[1]['logF_med_intensity_T'] > 0:
        new_row = pd.Series([row[1]['fragment'], row[1]['gene_names'],
                             row[1]['logF_med_intensity_T'], 'I'],
                            index=levy_df.columns)
        levy_df = levy_df.append(new_row, ignore_index=True)


levy_df['Subunit'] = "26S"
levy_df['full_sequence'] = None

# get unique protein ID's
annot['gene_upper'] = [x.upper() for x in annot['external_gene_name']]

split_gene_ids = [x.split(";") for x in levy_df['gene_names'].dropna()]
flattened_genes = [item for sublist in split_gene_ids for item in sublist]
unique_gene_ids = []
for item in flattened_genes:
    if item not in unique_gene_ids:
        unique_gene_ids.append(item)
unique_protein_ids = {}

for id in unique_gene_ids:
    prot_id = annot["ensembl_peptide_id"][annot["gene_upper"] == id]
    if prot_id.shape == (0,):
        continue
    else:
        unique_protein_ids[id] = prot_id.iloc[0]

# set up to store queried sequences, errors, and progress
sequence_dict = {}
error_index = []

# for each unique entry
print("Querying protein ID's")
for progress, entry in enumerate(unique_protein_ids):
    # try to query
    try:
        sequence_dict[entry] = retrieve_UniProt_seq(unique_protein_ids[entry])
    # if query fails, store index
    except:
        error_index.append(progress)
    # periodically output progress
    if progress % 100 == 0:
        print(round(progress/len(unique_protein_ids)*100, 3), "% completed")

assert len(error_index) == 0

levy_df = levy_df.dropna(subset=['gene_names'])
out_df = pd.DataFrame(columns=levy_df.columns)
out_df['frag_tracker'] = None

for i, entry in levy_df.iterrows():
    genes = entry['gene_names'].split(";")
    if len(genes) == 1:
        tmp_entry = pd.Series(entry)
        if tmp_entry['gene_names'] in sequence_dict.keys():
            tmp_entry['full_sequence'] = sequence_dict[tmp_entry['gene_names']]

        # append
        tmp_entry['frag_tracker'] = i
        out_df = out_df.append(tmp_entry, ignore_index=True)

    elif len(genes) > 1:
        tmp_df = pd.DataFrame(columns=out_df.columns)
        for gene in genes:
            tmp_entry = pd.Series(entry)
            tmp_entry['gene_names'] = gene
            if tmp_entry['gene_names'] in sequence_dict.keys():
                tmp_entry['full_sequence'] = sequence_dict[gene]
            tmp_df = tmp_df.append(tmp_entry, ignore_index=True)

        # append
        tmp_df['frag_tracker'] = i
        out_df = out_df.append(tmp_df, ignore_index=True)

    if i % 1000 == 0:
        print(i, "entry completed")

out_df = out_df.dropna(subset=["full_sequence"])

# TODO: add protein id's as full_accession_id column
# add columns needed in downstream processing
out_df['entry_source'] = "Levy_data"
out_df['start_pos'] = None
out_df['end_pos'] = None
out_df['origin_species'] = "human"
out_df['lit_reference'] = "https://doi.org/10.1038/nbt.4279"

out_df.to_csv(options.out, index=False)
