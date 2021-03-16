import os
import pandas as pd

in_dir = "/Users/weeder/PycharmProjects/pepsickle-paper/data/validation_data/" \
         "pt_data"
pt_files = []

for file in os.listdir(in_dir):
    if file.endswith(".csv"):
        pt_files.append(file)


ott_dat = pd.read_csv(in_dir + '/ott_peptides_with_context.csv')
ott_dat = ott_dat.drop_duplicates(['protein_context', 'response'])
ott_handle = in_dir + "/ott_peptides_with_context.fasta"
ott_val_fasta = open(ott_handle, "w")

for row in ott_dat.iterrows():
    index = row[0]
    entry = row[1]
    name = ">ott" + str(index) + "_" + str(entry['response'])
    ott_val_fasta.write(name)
    ott_val_fasta.write("\n")
    ott_val_fasta.write(entry['protein_context'])
    ott_val_fasta.write("\n")

ott_val_fasta.close()


mupexi_dat = pd.read_csv(in_dir + "/MuPeXI_peptides_with_context.csv")
mupexi_dat = mupexi_dat.drop_duplicates(['Immunogenic_Peptide',
                                         'protein_context'])
mupexi_handle = in_dir + '/mupexi_peptides_with_context.fasta'
mupexi_val_fasta = open(mupexi_handle, "w")

for row in mupexi_dat.iterrows():
    index = row[0]
    entry = row[1]
    name = ">mupexi" + str(index) + "_" + str(entry['Immunogenic_Peptide'])
    mupexi_val_fasta.write(name)
    mupexi_val_fasta.write("\n")
    mupexi_val_fasta.write(entry['protein_context'])
    mupexi_val_fasta.write("\n")

mupexi_val_fasta.close()


tesla_dat = pd.read_csv(in_dir + "/tesla_data_with_context.csv")
tesla_handle = in_dir + "/tesla_peptides_with_context.fasta"
tesla_val_fasta = open(tesla_handle, "w")

for row in tesla_dat.iterrows():
    index = row[0]
    entry = row[1]
    name = ">tesla" + str(index) + "_" + str(int(entry['response']))
    tesla_val_fasta.write(name)
    tesla_val_fasta.write("\n")
    tesla_val_fasta.write(entry['window'])
    tesla_val_fasta.write("\n")

tesla_val_fasta.close()
