#!usr/bin/env python3
"""

"""

import pandas as pd
from extraction_functions import *
from optparse import OptionParser


# define command line parameters
parser = OptionParser()
parser.add_option("-i", "--in_file", dest="in_file",
                  help="CSV of data from levy et al.")
parser.add_option("-o", "--out_dir", dest="out_dir",
                  help="output directory where processed CSV is exported")
parser.add_option("-a", "--annotation_file", dest="annotation_file",
                  help="CSV of ")

(options, args) = parser.parse_args()

# load in data
dat = pd.read_csv(options.in_file, low_memory=False, header=2, skiprows=[3],
                  usecols=range(0, 4))

annot = pd.read_csv(options.annotation_file, low_memory=False)

# fix column names
dat.columns = ['sequence', 'gene_names', 'logF_med_intensity_UT',
               'logF_med_intensity_T']

# remove decoy peptides
dat = dat[(dat['logF_med_intensity_T'].notna()) |
          (dat['logF_med_intensity_UT'].notna())]

out_df = pd.DataFrame(columns=['sequence', 'gene_names', 'intensity',
                               'Proteasome'])

for row in dat.iterrows():
    if row[1]['logF_med_intensity_UT'] > 0:
        new_row = pd.Series([row[1]['sequence'], row[1]['gene_names'],
                             row[1]['logF_med_intensity_UT'], 'C'],
                            index=out_df.columns)
        out_df = out_df.append(new_row, ignore_index=True)
    if row[1]['logF_med_intensity_T'] > 0:
        new_row = pd.Series([row[1]['sequence'], row[1]['gene_names'],
                             row[1]['logF_med_intensity_T'], 'I'],
                            index=out_df.columns)
        out_df = out_df.append(new_row, ignore_index=True)


out_df['Subunit'] = "26S"
out_df['full_sequence'] = None

# TODO: need to split gene lines with more than one entry... look up, and then determine if these are identical or not
# get unique protein ID's
annot['gene_upper'] = [x.upper() for x in annot['external_gene_name']]
unique_gene_ids = list(out_df['gene_names'].dropna().unique())
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
progress = 0

# for each unique entry
print("Querying protein ID's")
for entry in unique_protein_ids:
    # try to query
    try:
        sequence_dict[entry] = retrieve_UniProt_seq(entry)
    # if query fails, store index
    except:
        error_index.append(progress)
    progress += 1
    # periodically output progress
    if progress % 10 == 0:
        print(round(progress/len(unique_protein_ids)*100, 3), "% completed")


# attempt to repair null queries
progress = 0
print("Attempting to repair failed queries")
for e in error_index:
    # pull id of failed query
    tmp_id = unique_protein_ids[e]
    # ignore those with no real id
    if tmp_id != "not applicable":
        try:
            # re-query with broader parameters
            query = compile_UniProt_url(tmp_id, include_experimental=True)
            buffer = extract_UniProt_table(query)
            new_id = buffer["Entry"][0]
            sequence_dict[unique_protein_ids[e]] = retrieve_UniProt_seq(new_id)
        except IndexError:
            # if empty results table, pass
            pass
    progress += 1
    # periodically output progress
    if progress % 10 == 0:
        print(round(progress/len(error_index)*100, 3), "% completed")


# for each entry in the results table
for e in range(len(out_df)):
    # pull protein ID
    prot_id = str(out_df.at[e, 'gene_names'])
    # if sequence was found for given ID, store full sequence
    if prot_id in sequence_dict.keys():
        out_df.at[e, 'full_sequence'] = sequence_dict[prot_id]

# drop rows where no sequence was retrieved
out_df.dropna(subset=['full_sequence'], inplace=True)

# add columns needed in downstream processing
out_df['entry_source'] = "Levy_data"
out_df['start_pos'] = None
out_df['end_pos'] = None


