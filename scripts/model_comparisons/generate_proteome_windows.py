from Bio import SeqIO
from random import sample
from sequence_featurization_tools import *
import pandas as pd

fasta_file = "/Users/weeder/PycharmProjects/pepsickle-paper/data/raw/" \
          "human_proteome.fasta"

proteins = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
unique_windows = set()

for id in proteins:
    tmp_peptide = str(proteins[id].seq)
    peptide_windows = create_windows_from_protein(tmp_peptide, upstream=10,
                                                  downstream=10)
    unique_windows.update(peptide_windows)

# sample_windows = sample(unique_windows, round(len(unique_windows)*.25))
sample_windows = unique_windows

proteome_window_features = generate_feature_array(sample_windows)
print(proteome_window_features.shape)

m, n, r = proteome_window_features[:, :, 22:].shape
tmp_dat = proteome_window_features[:, :, 22:].reshape(-1, n*r)
print(tmp_dat.shape)
dat = pd.DataFrame(tmp_dat)
dat['id'] = "human_proteome"


dat.to_csv("/Users/weeder/PycharmProjects/pepsickle-paper/data/validation_data/"
           "output/plots/whole_proteome_frag_features_21aa.csv")
