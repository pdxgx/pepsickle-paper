from Bio import SeqIO
from collections import defaultdict

# Create gene to protein sequences and gene/protein sequences to protein ID dictionaries
gene_to_proteins_dict = defaultdict(set)
protein_to_id_dict = {}
for record in SeqIO.parse('/Users/weeder/PycharmProjects/pepsickle-paper/data/validation_data/pt_data/tmp/gencode.v36.pc_translations.fa', 'fasta'):
	gene = str(record.id).split('|')[6]
	protein = str(record.seq)
	protein_id = str(record.id).split('|')[0].strip('>')
	gene_to_proteins_dict[gene].add(protein)
	protein_to_id_dict[(gene, protein)] = protein_id


# Process peptides
with open('/Users/weeder/PycharmProjects/pepsickle-paper/data/validation_data/pt_data/ott_peptides.csv') as f:
	with open('/Users/weeder/PycharmProjects/pepsickle-paper/data/validation_data/pt_data/ott_peptides_with_context.csv', 'w') as o:
		print(','.join(['epitope', 'protein_context', 'aa_change', 'protein_id', 'response']), file=o)
		for i in range(4):
			f.readline()
		for line in f:
			tokens = line.strip().split(',')
			normal_pep = tokens[8]
			# Check whether valid peptide:
			if normal_pep not in ['', 'NA'] and tokens[14] in ['0', '1']:
				mutant_pep = tokens[6]
				#print(normal_pep, mutant_pep)
				# Get aa change positional information
				aa_change = tokens[3]
				position = int(aa_change.split('.')[1][1:-1]) - 1
				change_pos = [i for i in range(len(mutant_pep)) if mutant_pep[i] != normal_pep[i]][0]
				# Get 0-based start position of peptide in protein:
				start_pos = position - change_pos
				end_pos = start_pos + len(normal_pep) - 1
				#print(aa_change, change_pos, len(mutant_pep), start_pos, end_pos)
				# Get protein context
				gene = tokens[2].rstrip('a').rstrip('b').rstrip('c')
				found_match = False
				for protein in gene_to_proteins_dict[gene]:
					# Check if normal protein is the appropriate position
					if protein[start_pos:start_pos+len(normal_pep)] == normal_pep:
						# Get sequence context
						left_context = protein[max(0, start_pos-16):start_pos]
						if len(left_context) < 16:
							leftover = 16 - len(left_context)
							left_context = ''.join([''.join(['*' for i in range(leftover)]), left_context])
						right_context = protein[end_pos+1:end_pos+9]
						if len(right_context) < 8:
							leftover = 8 - len(right_context)
							right_context = ''.join([right_context, ''.join(['*' for i in range(leftover)])])
						full_sequence = ''.join([left_context, mutant_pep, right_context])
						#print(left_context, right_context, len(protein))
						assert len(full_sequence) == len(mutant_pep) + 24
						# Get annotation info and write output
						protein_id = protein_to_id_dict[(gene, protein)]
						found_match = True
						print(','.join([mutant_pep, full_sequence, aa_change, protein_id, tokens[14]]), file=o)
				if not found_match:
					print(line.strip())
					print(gene_to_proteins_dict[gene])
					print()
