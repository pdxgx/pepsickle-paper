from Bio import SeqIO
from Bio.Blast import NCBIWWW
import subprocess
import pandas as pd
import os
import argparse
import subprocess
import os
import datetime
import pickle
import re
from Bio.SubsMat import MatrixInfo


def process_blast(blast_results, type, matrix, dict_dir):
    ''' Processes results of blastp to obtain dictionary of best results
        Keys are neoepitope sequences
        Values are lists of associated blast data: E value, sequence of match, raw protein similarity score
        For human blast, transcript and gene of match peptide are also stored
        For bacterial or viral blast, species of match peptide is also stored
        blast_results: path to file containing results of blastp
        type: blast type - human, bacterial, or viral (string)
        matrix: scoring matrix to use for comparing sequences (matrix)
        Return value: dictionary
    '''
    # Set peptide dictionary based on peptide type
    if type == "human":
        dict = pickle.load(open(dict_dir + "humanDict.pickle", "rb"))
    elif type == "bacterial":
        dict = pickle.load(open(dict_dir + "bacterialDict.pickle", "rb"))
    elif type == "viral":
        dict = pickle.load(open(dict_dir + "viralDict.pickle", "rb"))
    # Parse blast output line-by-line
    blast_dict = {}
    with open(blast_results, "r") as fh:
        for line in fh:
            # Obtain relevant data - epitope sequence, alignment length, E value, matching sequence, peptide name
            line = line.strip("\n").split("\t")
            epitope = line[0].split("=")[1]
            length = int(line[2])
            eval = float(line[6])
            match_seq = line[5]
            match_pep = line[1].replace("ref", "").replace("|", "")
            # If working with human peptides, obtain match transcript and match gene from dictionary
            if type == "human":
                match_transcript = dict[match_pep][1]
                match_gene = dict[match_pep][0]
            # If working with bacterial or viral peptides, obtain the match species
            else:
                match_species = dict[match_pep]
            # Check for presence of invalid characters in match seq
            invalids = ["B", "J", "O", "U", "X", "Z", "*"]
            invalid_matches = []
            for char in invalids:
                if char in match_seq:
                    invalid_matches.append(char)
            # Check f epitope is not already in dictionary, the alignment is the right length, and there are no invalid characters
            if epitope not in blast_dict and length == len(epitope) and invalid_matches == []:
                # Find BLOSUM score between neoepitope and peptide match, then store data
                match_ps = score_pairwise(epitope, match_seq, matrix)
                if type == "human":
                    blast_dict[epitope] = [eval, match_transcript, match_gene, match_seq, match_ps]
                else:
                    blast_dict[epitope] = [eval, match_species, match_seq, match_ps]
            # If epitope is in dictionary, but E value for this entry is better, replace data
            elif epitope in blast_dict and length == len(epitope) and eval < blast_dict[epitope][
                0] and invalid_matches == []:
                # Find BLOSUM score between neoepitope and peptide match, then store data
                match_ps = score_pairwise(epitope, match_seq, matrix)
                if type == "human":
                    blast_dict[epitope] = [eval, match_transcript, match_gene, match_seq, match_ps]
                else:
                    blast_dict[epitope] = [eval, match_species, match_seq, match_ps]
            # If epitope is in dictionary and E value for this entry is equivalent, compare further
            # If the match sequence is the same as previous entry, store this entry's data too
            # If match sequence is a better match to neoepitope, replace data
            elif epitope in blast_dict and length == len(epitope) and eval == blast_dict[epitope][
                0] and invalid_matches == []:
                if type == "human":
                    if match_seq == blast_dict[epitope][3] and match_transcript not in blast_dict[epitope][1] \
                            and match_gene not in blast_dict[epitope][2]:
                        blast_dict[epitope][1] = blast_dict[epitope][1] + "," + match_transcript
                        blast_dict[epitope][2] = blast_dict[epitope][2] + "," + match_gene
                    else:
                        match_ps = score_pairwise(epitope, match_seq, matrix)
                        if match_seq == epitope or match_ps > blast_dict[epitope][4]:
                            blast_dict[epitope] = [eval, match_transcript, match_gene, match_seq, match_ps]
                else:
                    if match_seq == blast_dict[epitope][2] and match_species not in blast_dict[epitope][1]:
                        blast_dict[epitope][1] = blast_dict[epitope][1] + "," + match_species
                    else:
                        match_ps = score_pairwise(epitope, match_seq, blosum)
                        if match_seq == epitope or match_ps > blast_dict[epitope][3]:
                            blast_dict[epitope] = [eval, match_species, match_seq, match_ps]
    return blast_dict


def run_blast(fasta, db, blastp, outputdir, name, type):
    ''' Runs blastp
        fasta: path fasta containing sequences to blast against database
        db: path to database to blast against
        outputdir: path to directory in which to write blast output
        name: sample name to distinguish file (string)
        type: blast type - human, bacterial, or viral (string)
        No return value
    '''
    outfile = outputdir + "/" + name + "." + type + ".blast.out"
    # If blast output does not already exist, run blast to obtain it
    if os.path.isfile(outfile) == False:
        subprocess.call(
            [blastp, "-outfmt", "6 qseqid sseqid length nident sstart send", "-db", db, "-query", fasta, "-matrix",
             "BLOSUM62", "-evalue", "200000", "-ungapped", "-comp_based_stats", "F", "-out", outfile])
    # If blast output already exits, skip running blast
    else:
        print
        "Blast file for " + type + " peptides already exists - skipping"


def make_epitope_fasta(epitope_file, outputdir, name, fasta):
    ''' Produces fasta file containing all neoepitope sequences for a sample
        epitope_file: path to parsed output file from pVAC-Seq (or alternative program)
        outputdir: path to directory in which to write fasta
        name: sample name to distinguish file (string)
        fasta: path to output fasta file
        No return value
    '''
    # Obtain all unique epitopes and store as a list
    epitope_list = []
    with open(epitope_file, "r") as fh:
        for line in fh:
            line = line.split("\t")
            epitope = line[1]
            if epitope not in epitope_list:
                epitope_list.append(epitope)
    # Write unique epitopes from list to fasta file
    with open(fasta, "w") as fh:
        for epitope in epitope_list:
            fh.write("> seq=" + epitope + "\n")
            fh.write(epitope + "\n")


in_dir = '/Users/weeder/PycharmProjects/pepsickle-paper/data/validation_data/pt_data'
tesla_dat = pd.read_excel(in_dir + "/1-s2.0-S0092867420311569-mmc4.xlsx")
tesla_dat = tesla_dat.dropna(axis=0, subset=['MUTATION_POSITION'])

tesla_handle = in_dir + "/tesla_epitopes.fasta"
tesla_out = open(tesla_handle, "w")

for row in tesla_dat.iterrows():
    i, entry = row
    index = row[0]
    entry = row[1]
    mut_ep = str(entry['ALT_EPI_SEQ'])
    mut_pos = int(entry['MUTATION_POSITION'])
    name = ">" + str(mut_ep) + "_" + str(mut_pos) + "_" + str(entry['VALIDATED'])
    base_ep = mut_ep[:mut_pos-1] + "X" + mut_ep[mut_pos:]

    tesla_out.write(name)
    tesla_out.write("\n")
    tesla_out.write(entry['ALT_EPI_SEQ'])
    tesla_out.write("\n")

tesla_out.close()

humanDB = "/Users/weeder/PycharmProjects/pepsickle-paper/data/validation_data/pt_data/blast_dbs/hg38_peptide_db"
run_blast(fasta=in_dir+"/tesla_epitopes.fasta", db=humanDB, outputdir=in_dir,
          name="tesla_epitopes", type="human", blastp="/usr/local/ncbi/blast/bin/blastp")

blast_in = in_dir + "/tesla_epitopes.human.blast.out"
blast_dat = pd.read_csv(blast_in, sep="\t", header=None)
unique_query = blast_dat[0].unique()

source_dict = {}
query_errors = []
mismatch_errors = []
position_errors = []
counter = 0
for query in unique_query:
    print(round(counter/len(unique_query), 3))
    counter += 1

    sub_blast = blast_dat[blast_dat[0] == query]
    orig_ep = query.split("_")[0]
    sub_blast = sub_blast[((sub_blast[2] == len(orig_ep)) &
                           (sub_blast[3] == len(orig_ep)-1)) |
                          ((sub_blast[2] == len(orig_ep)-1) &
                           (sub_blast[3] == len(orig_ep)-1))]

    if len(sub_blast) == 0:
        query_errors.append(query)
    else:
        for row in sub_blast.iterrows():
            i, entry = row
            p_id = entry[1].split(".")[0]
            prot_source = retrieve_UniProt_seq(p_id)
            search_string = orig_ep[:int(query.split("_")[1])-1] + "." + orig_ep[int(query.split("_")[1]):]
            start = entry[4] - 1
            end = entry[5]

            if len(prot_source) > 0:
                try:
                    search_match = re.search(search_string, prot_source).span()
                except AttributeError:
                    if query not in source_dict.keys() and query not in mismatch_errors:
                        mismatch_errors.append(query)
                    continue

                if int(query.split("_")[1]) == 1 and len(query.split("_")[0])-1 == int(entry[2]):
                    start = int(entry[4]) - 1
                if int(query.split("_")[1]) == len(query.split("_")[0]) and len(query.split("_")[0])-1 == int(entry[2]):
                    end = int(entry[5]) + 1

                if prot_source[start:end] == prot_source[search_match[0]:search_match[1]]:
                    if query not in source_dict.keys():
                        source_dict[query] = []
                    if prot_source not in source_dict[query]:
                        source_dict[query].append(prot_source)
                else:
                    if query not in position_errors and query not in source_dict.keys():
                        position_errors.append(query)
            else:
                if query not in query_errors:
                    query_errors.append(query)


successful_entries = []
for query in unique_query:
    if query in source_dict.keys():
        successful_entries.append(query)
    if query not in query_errors and query not in mismatch_errors and query not in position_errors:
        if query not in source_dict.keys():
            print(query)


out_df = pd.DataFrame()
ambiguous_windows = []
for entry in successful_entries:
    orig_ep = entry.split("_")[0]
    if entry.split("_")[2] == "False":
        response = int(0)
    else:
        response = int(1)

    search_string = orig_ep[:int(entry.split("_")[1])-1] + "." + orig_ep[int(entry.split("_")[1]):]
    possible_sources = source_dict[entry]

    if len(possible_sources) == 1:
        source = possible_sources[0]
        match = re.search(search_string, source)
        start = match.span()[0]-16
        end = match.span()[1] + 8

        if start < 0:
            start_buffer = "*" * abs(start)
            start = 0
        else:
            start_buffer = ""

        if end > len(source):
            end_buffer = "*" * (end - len(source))
            end = len(source)
        else:
            end_buffer = ""

        window = start_buffer + source[start:match.span()[0]] + \
                 orig_ep + source[match.span()[1]:end] + end_buffer
        row = pd.Series([orig_ep, window, response])
        out_df = out_df.append(row, ignore_index=True)

    else:
        possible_windows = []
        for s in possible_sources:
            source = s
            match = re.search(search_string, source)
            start = match.span()[0]-16
            end = match.span()[1] + 8

            if start < 0:
                start_buffer = "*" * abs(start)
                start = 0
            else:
                start_buffer = ""

            if end > len(source):
                end_buffer = "*" * (end - len(source))
                end = len(source)
            else:
                end_buffer = ""

            p_window = start_buffer + source[start:match.span()[0]] + \
                     orig_ep + source[match.span()[1]:end] + \
                     end_buffer

            if str(p_window) not in possible_windows:
                possible_windows.append(p_window)

        if len(possible_windows) == 1:
            window = possible_windows[0]
            row = pd.Series([orig_ep, window, response])
            out_df = out_df.append(row, ignore_index=True)
        else:
            print(possible_windows)
            ambiguous_windows.append(entry)
        # print(len(possible_windows))

out_df.columns = ['epitope', 'window', 'response']
out_df.to_csv(in_dir + "/tesla_data_with_context.csv", index=False)
