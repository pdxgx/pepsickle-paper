#!/usr/bin/env python3
"""
[].py

For issues contact Ben Weeder (weeder@ohsu.edu)

"""

import pandas as pd
import re
from extraction_functions import *
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--in-file",
                  help="pickled dictionary of cleavage windows")
parser.add_option("-o", "--out",
                  help="output file where results will be exported")

(options, args) = parser.parse_args()

# load data
dat = pd.read_csv(options.in_file, low_memory=False, skiprows=0, header=0)

# extract unique id's
# unique_gene_ids = dat['Gene_ID'].unique()
unique_source_sequences = {}
failed_ids = []

unique_tx_ids = []
for tx_id in dat['Transcript_ID']:
    tmp_ids = tx_id.split(",")
    for new_id in tmp_ids:
        if new_id not in unique_tx_ids:
            unique_tx_ids.append(new_id)

for id in unique_tx_ids:
    # print(id)
    id_query = compile_UniProt_url(id)
    try:
        prot_id = extract_UniProt_table(id_query)['Entry'][0]
        unique_source_sequences[id] = retrieve_UniProt_seq(prot_id)
    except:
        failed_ids.append(id)

print("Failed queries:", len(failed_ids), "attempting to repair...")
for id in failed_ids:
    id_query = compile_UniProt_url(id, include_experimental=True)
    try:
        prot_id = extract_UniProt_table(id_query)['Entry'][0]
        unique_source_sequences[id] = retrieve_UniProt_seq(prot_id)
    except:
        print(id, "Failed")


# to revert to source ep
test_row = dat.iloc[0]
test_ep = list(test_row["Mut_peptide"])
test_ep[test_row['peptide_position']-1] = test_row['Amino_Acid_Change'][0]
base_ep = "".join(test_ep)

base_epitopes = []
for row in dat.iterrows():
    mut_ep = list(row[1]['Mut_peptide'])
    mut_ep[row[1]['peptide_position']-1] = row[1]['Amino_Acid_Change'][0]
    base_ep = "".join(mut_ep)
    base_epitopes.append(base_ep)

dat['base_peptide'] = base_epitopes

out_df = pd.DataFrame()
for row in dat.iterrows():
    entry = row[1]
    tx_list = entry['Transcript_ID'].split(",")

    source_seqs = []
    for tx in tx_list:
        try:
            source_seqs.append(unique_source_sequences[tx])
        except KeyError:
            pass

    tmp_windows = []
    for seq in source_seqs:
        matches = re.finditer(entry['base_peptide'], seq)
        for match in matches:
            start = match.span()[0]-16
            end = match.span()[1] + 8

            if start < 0:
                start_buffer = "*" * abs(start)
                start = 0
            else:
                start_buffer = ""

            if end > len(seq):
                end_buffer = "*" * (end - len(seq))
                end = len(seq)
            else:
                end_buffer = ""
            window = start_buffer + seq[start:match.span()[0]] + \
                     entry['Mut_peptide'] + seq[match.span()[1]:end] + \
                     end_buffer

            if window not in tmp_windows:
                tmp_windows.append(window)

    static_vals = [entry["Gene_ID"], entry['Mut_peptide'],
                   entry['Amino_Acid_Change'],
                   entry['priority_Score'], entry['Pateint'],
                   entry['Ordered_Peptide'], entry['Immunogenic_Peptide']]

    for w in tmp_windows:
        out_row = pd.Series(static_vals+[w])
        out_df = out_df.append(out_row, ignore_index=True)

out_df.columns = ["Gene_ID", 'Mut_peptide', 'Amino_Acid_Change',
                  'priority_Score', 'Pateint', 'Ordered_Peptide',
                  'Immunogenic_Peptide', 'protein_context']

print("Exporting Final Table: ")
out_df.to_csv(options.out, index=False)
