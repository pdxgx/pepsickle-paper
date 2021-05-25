#!bin/bash

#### SETUP
## set working directory to base dir of project
cd /Users/weeder/PycharmProjects/pepsickle-paper
## set temp environmental vars for mysql use
export MYSQL_USER=[USER]
export MYSQL_PWD=[PASSWORD]

#### pull database data
## ANTIJEN
## NOTE: issue with antijen database... database is no longer available. See static data extract for use in downstream processes
python3 ./scripts/database_pulls/extract_AntiJen_data.py -o ./data/raw/database_pulls/
## SYFPEITHI
python3 ./scripts/database_pulls/extract_SYF_data.py -a ./data/raw/raw_mammal_allele_series.csv -o ./data/raw/database_pulls/

#### extract static data
## IEDB
wget -O ./data/raw/database_pulls/iedb_public.sql.gz https://www.iedb.org/downloader.php?file_name=doc/iedb_public.sql.gz
gzip -d ./data/raw/database_pulls/iedb_public.sql.gz
mysql --password=$MYSQL_PWD -u $MYSQL_USER < ./scripts/static_dataset_extractions/iedb_mysql_query.sql
## extract BC data
python3 ./scripts/static_dataset_extractions/extract_breast_cancer_data.py -i ./data/raw/breast_cancer_data/ -o ./data/processed/
## extract digestion data (in-vitro study data)
python3 ./scripts/static_dataset_extractions/extract_digestion_data.py -i ./data/raw/digestion_map_files -o ./data/processed/

#### process extracted and static data
## database data
# AntiJen
python3 ./scripts/dataset_processing/get_AntiJen_source_sequences.py -i ./data/raw/database_pulls/AntiJen_Tcell_epitopes.csv -o ./data/processed/
# SYFPEITHI
python3 ./scripts/dataset_processing/get_SYFPEITHI_source_sequences.py -i ./data/raw/database_pulls/SYFPEITHI_epitopes.csv -o ./data/processed/
# IEDB
python3 ./scripts/dataset_processing/iedb_mysql_processing.py -u $MYSQL_USER -p $MYSQL_PWD -o ./data/processed/
# Breast cancer data
python3 ./scripts/dataset_processing/get_breast_cancer_source_sequences.py -i ./data/processed/breast_cancer_epitopes.csv -o ./data/processed/
# Winter et al
python3 ./scripts/dataset_processing/prep_winter_data.py -i ./data/raw/Winter_et_al_results.csv -o ./data/processed/

#### Merge datasets together
python3 ./scripts/merging_and_filtering/merge_datasets.py -i ./data/processed -o ./data/merged
## NOTE: run if only human data is desired
# python3 ./scripts/merging_and_filtering/merge_datasets.py -i ./data/processed -o ./data/merged --human-only

#### Verify indices in full file
python3 ./scripts/merging_and_filtering/epitope_index_check.py -i ./data/merged/merged_data_all_mammal.csv -o ./data/merged/merged_data_all_mammal_verified.csv
# if only human data is desired
# python3 ./scripts/merging_and_filtering/epitope_index_check.py -i ./data/merged/merged_data_human_only.csv -o ./data/merged/merged_data_human_only_verified.csv

#### generate the negative fragment examples based on annotated positives
## NOTE: final window sized are shown, but -w option allows for variable window sizes
## for 20S
python3 ./scripts/merging_and_filtering/negative_set_generation.py --full-negative-set -w 3,3 -i ./data/merged/merged_data_all_mammal_verified.csv -o ./data/training_sets/all_mammal_windows_7aa.pickle
## for epitopes
python3 ./scripts/merging_and_filtering/negative_set_generation.py --full-negative-set -w 8,8 -i ./data/merged/merged_data_all_mammal_verified.csv -o ./data/training_sets/all_mammal_windows_17aa.pickle

#### run model training
## NOTE: uncommented script was used to train the most performant model
## 20S digestion models
# NN model
# python3 ./scripts/modeling/digestion_map_based_ensemble_net.py -w 7 -i ./data/training_sets/all_mammal_windows_7aa.pickle -o ./data/model_weights/window_7aa
# ML model
## NOTE: when using chemical features, values are normalized
python3 ./scripts/modeling/proteasome_GBM_model.py -m proteasome -f chemistry -n -i ./data/training_sets/all_mammal_windows_7aa.pickle -o ./data/model_weights/ML_weights/in-vitro_7aa/chem
# python3 ./scripts/modeling/proteasome_GBM_model.py -m proteasome -f identity -i ./data/training_sets/all_mammal_windows_7aa.pickle -o ./data/model_weights/ML_weights/in-vitro_7aa/identity

## epitope models
python3 ./scripts/modeling/epitope_based_ensemble_net.py -w 17 -i ./data/training_sets/all_mammal_windows_17aa.pickle -o ./data/model_weights/DL_weights/window_17aa
# ML model
# python3 ./scripts/modeling/proteasome_GBM_model.py -m epitope -f chemistry -n -i ./data/training_sets/all_mammal_windows_17aa.pickle -o ./data/model_weights/ML_weights/epitope_7aa/chemistry
# python3 ./scripts/modeling/proteasome_GBM_model.py -m epitope -f identity  -i ./data/training_sets/all_mammal_windows_17aa.pickle -o ./data/model_weights/ML_weights/epitope_17aa/identity

#### compile model weights into single file
## NOTE: if you'd like to use models other than the primary models for end-point analyses, move all generated weights to the same folder before running this command
python3 ./scripts/modeling/generate_model_dict.py -i ./data/model_weights/DL_weights/window_17aa -o ./data/model_weights/DL_weights/window_17aa

#### Prepare validation set
## compile fragments for left out 20S digestion validation_prep data
python3 ./scripts/static_dataset_extractions/extract_digestion_data.py -i ./data/validation_data/digestion_data/raw/ -o ./data/validation_data/digestion_data/
python3 ./scripts/validation_prep/prep_20S_digestion_val_data.py -i ./data/validation_data/digestion_data/compiled_digestion_df.csv -o ./data/validation_data/digestion_data/
## NOTE: --full-negative not used to give balance initial set
python3 ./scripts/merging_and_filtering/epitope_index_check.py -i ./data/validation_data/digestion_data/20S_digestion_val_columns_remapped.csv -o ./data/validation_data/digestion_data/20S_digestion_val_data_indices_verified.csv

## compile epitope validation_prep data
python3 ./scripts/validation_prep/prep_epitope_validation_data.py -i ./data/validation_data/epitope_data/41467_2016_BFncomms13404_MOESM1318_ESM.csv -o ./data/validation_data/epitope_data
python3 ./scripts/merging_and_filtering/epitope_index_check.py -i ./data/validation_data/epitope_data/validation_epitopes_w_source.csv -o ./data/validation_data/epitope_data/validation_epitopes_verified.csv

#### convert entries to cleavage windows
## NOTE: largest window size is selected w/ downsampling for unique internal windows downstream
python3 ./scripts/merging_and_filtering/negative_set_generation.py -w 9,9 -i ./data/validation_data/digestion_data/20S_digestion_val_data_indices_verified.csv -o ./data/validation_data/validation_sets_pre-filter/20S_val_windows_21aa_paired.pickle
python3 ./scripts/merging_and_filtering/negative_set_generation.py -w 9,9 -i ./data/validation_data/epitope_data/validation_epitopes_verified.csv -o ./data/validation_data/validation_sets_pre-filter/epitope_val_windows_21aa_paired.pickle


## filter out entries seen in training data
python3 ./scripts/validation_prep/filter_val_data.py -w 21 --internal-filter-size 7 -i ./data/validation_data/validation_sets_pre-filter -t ./data/training_sets -o ./data/validation_data/completed_validation_sets

#### Install pepsickle to handle application of the models to validation data. The full instructions can be found at [URL]
## NOTE: to use the weights directly here, simply replace the weight files (model.joblib and trained_model_dict.pickle) that come with pepsickle by default.

#### Run predictions on validation data
pepsickle -v -f ./data/validation_data/completed_validation_sets/window_fasta_files/epitope_val_data_21aa.fasta -o ./data/validation_data/validation_results/pepsickle/pepsickle_epitope_val_results.txt
pepsickle -m in-vitro -p C -f ./data/validation_data/completed_validation_sets/window_fasta_files/20S_digestion_constitutive_validation_data_21aa.fasta -o ./data/validation_data/validation_results/pepsickle/pepsickle_constitutive_20S_val_results.txt
pepsickle -m in-vitro -p I -f ./data/validation_data/completed_validation_sets/window_fasta_files/20S_digestion_immuno_validation_data_21aa.fasta -o ./data/validation_data/validation_results/pepsickle/pepsickle_immuno_20S_val_results.txt

## summarize pepsickle results
python3 ./scripts/model_comparisons/process_pepsickle_results.py -i ./data/validation_data/validation_results/pepsickle/pepsickle_epitope_val_results.txt -o ./data/validation_data/validation_results/pepsickle/pepsickle_epitope_val_summary.txt
python3 ./scripts/model_comparisons/process_pepsickle_results.py -i ./data/validation_data/validation_results/pepsickle/pepsickle_constitutive_20S_val_results.txt -o ./data/validation_data/validation_results/pepsickle/pepsickle_constitutive_20S_val_summary.txt
python3 ./scripts/model_comparisons/process_pepsickle_results.py -i ./data/validation_data/validation_results/pepsickle/pepsickle_immuno_20S_val_results.txt -o ./data/validation_data/validation_results/pepsickle/pepsickle_immuno_20S_val_summary.txt


#### Compare to other cleavage tools
## after running PCPS on val fastas... summarize results
python3 ./scripts/model_comparisons/process_pcps_results.py -i ./data/validation_data/validation_results/PCPS/PCPS_epitope_preds.csv -p 21 -o ./data/validation_data/validation_results/PCPS/PCPS_epitope_pred_summary.csv
python3 ./scripts/model_comparisons/process_pcps_results.py -i ./data/validation_data/validation_results/PCPS/PCPS_20S_constit_preds.csv -p 21 -o ./data/validation_data/validation_results/PCPS/PCPS_20S_constit_pred_summary.csv
python3 ./scripts/model_comparisons/process_pcps_results.py -i ./data/validation_data/validation_results/PCPS/PCPS_20S_immuno_preds.csv -p 21 -o ./data/validation_data/validation_results/PCPS/PCPS_20S_immuno_summary.csv

## after running NetChop on val fastas... summarize
python3 ./scripts/model_comparisons/parse_netchop_file.py -i ./data/validation_data/validation_results/netchop/netchop_epitope_21aa_val_cleavage_preds.txt -o ./data/validation_data/validation_results/netchop/parsed_netchop_epitope_val_21aa_cleavage_preds.txt
python3 ./scripts/model_comparisons/parse_netchop_file.py -i ./data/validation_data/validation_results/netchop//netchop_20S_constitutive_digestion_val_21aa_cleavage_preds.txt -o ./data/validation_data/validation_results/netchop/parsed_netchop_20S_constitutive_digestion_val_21aa_cleavage_preds.txt
python3 ./scripts/model_comparisons/parse_netchop_file.py -i ./data/validation_data/validation_results/netchop/netchop_20S_immuno_digestion_val_21aa_cleavage_preds.txt -o ./data/validation_data/validation_results/netchop/parsed_netchop_20S_immuno_digestion_val_21aa_cleavage_preds.txt


## cross performance
pepsickle -v -f ./data/validation_data/completed_validation_sets/window_fasta_files/20S_digestion_constitutive_validation_data_21aa.fasta -o ./data/validation_data/validation_results/pepsickle/pepsickle_epitope_cross_constit_val_results.txt
pepsickle -v -f ./data/validation_data/completed_validation_sets/window_fasta_files/20S_digestion_immuno_validation_data_21aa.fasta -o ./data/validation_data/validation_results/pepsickle/pepsickle_epitope_cross_immuno_val_results.txt
pepsickle -m in-vitro -p C -f ./data/validation_data/completed_validation_sets/window_fasta_files/epitope_val_data_21aa.fasta -o ./data/validation_data/validation_results/pepsickle/pepsickle_constitutive_20S_cross_epitope_val_results.txt
pepsickle -m in-vitro -p I -f ./data/validation_data/completed_validation_sets/window_fasta_files/epitope_val_data_21aa.fasta -o ./data/validation_data/validation_results/pepsickle/pepsickle_immuno_20S_cross_epitope_val_results.txt
## parse cross performance
python3 ./scripts/model_comparisons/process_pepsickle_results.py -i ./data/validation_data/validation_results/pepsickle/pepsickle_constitutive_20S_cross_val_results.txt -o ./data/validation_data/validation_results/pepsickle/pepsickle_constitutive_20S_cross_val_summary.txt
python3 ./scripts/model_comparisons/process_pepsickle_results.py -i ./data/validation_data/validation_results/pepsickle/pepsickle_immuno_20S_cross_val_results.txt -o ./data/validation_data/validation_results/pepsickle/pepsickle_immuno_20S_cross_val_summary.txt

python3 ./scripts/model_comparisons/process_pepsickle_results.py -i ./data/validation_data/validation_results/pepsickle/pepsickle_epitope_cross_constit_val_results.txt -o ./data/validation_data/validation_results/pepsickle/pepsickle_epitope_cross_constit_val_summary.txt
python3 ./scripts/model_comparisons/process_pepsickle_results.py -i ./data/validation_data/validation_results/pepsickle/pepsickle_epitope_cross_immuno_val_results.txt -o ./data/validation_data/validation_results/pepsickle/pepsickle_epitope_cross_immuno_val_summary.txt

# cross performance iv <-> constit
pepsickle -m in-vitro -p C -f ./data/validation_data/completed_validation_sets/window_fasta_files/20S_digestion_immuno_validation_data_21aa.fasta -o ./data/validation_data/validation_results/pepsickle/pepsickle_constitutive_20S_cross_val_immuno_results.txt
pepsickle -m in-vitro -p I -f ./data/validation_data/completed_validation_sets/window_fasta_files/20S_digestion_constitutive_validation_data_21aa.fasta -o ./data/validation_data/validation_results/pepsickle/pepsickle_immuno_20S_cross_val_constit_results.txt
## parse cross performance
python3 ./scripts/model_comparisons/process_pepsickle_results.py -i ./data/validation_data/validation_results/pepsickle/pepsickle_constitutive_20S_cross_val_immuno_results.txt -o ./data/validation_data/validation_results/pepsickle/pepsickle_constitutive_20S_cross_val_immuno_summary.txt
python3 ./scripts/model_comparisons/process_pepsickle_results.py -i ./data/validation_data/validation_results/pepsickle/pepsickle_immuno_20S_cross_val_constit_results.txt -o ./data/validation_data/validation_results/pepsickle/pepsickle_immuno_20S_cross_val_constit_summary.txt


## run on pt data (epitope model)
pepsickle -f ./data/validation_data/pt_data/tesla_peptides_with_context.fasta -o ./data/validation_data/pt_data/tesla_cleavage_preds.txt
pepsickle -f ./data/validation_data/pt_data/ott_peptides_with_context.fasta -o ./data/validation_data/pt_data/ott_cleavage_preds.txt
pepsickle -f ./data/validation_data/pt_data/mupexi_peptides_with_context.fasta -o ./data/validation_data/pt_data/mupexi_cleavage_preds.txt
## run on pt data (iv model)
pepsickle -m in-vitro -p I -f ./data/validation_data/pt_data/tesla_peptides_with_context.fasta -o ./data/validation_data/pt_data/tesla_cleavage_preds_iv.txt
pepsickle -m in-vitro -p I -f ./data/validation_data/pt_data/ott_peptides_with_context.fasta -o ./data/validation_data/pt_data/ott_cleavage_preds_iv.txt
pepsickle -m in-vitro -p I -f ./data/validation_data/pt_data/mupexi_peptides_with_context.fasta -o ./data/validation_data/pt_data/mupexi_cleavage_preds_iv.txt

python3 ./scripts/model_comparisons/pepsickle_pt_data_processing.py --in-vivo ./data/validation_data/pt_data/tesla_cleavage_preds.txt --in-vitro ./data/validation_data/pt_data/tesla_cleavage_preds_iv.txt -o ./data/validation_data/pt_data/pepsickle_pt_summaries/tesla_pt_summary.csv
python3 ./scripts/model_comparisons/pepsickle_pt_data_processing.py --in-vivo ./data/validation_data/pt_data/ott_cleavage_preds.txt --in-vitro ./data/validation_data/pt_data/ott_cleavage_preds_iv.txt -o ./data/validation_data/pt_data/pepsickle_pt_summaries/ott_pt_summary.csv
python3 ./scripts/model_comparisons/pepsickle_pt_data_processing.py --in-vivo ./data/validation_data/pt_data/mupexi_cleavage_preds.txt --in-vitro ./data/validation_data/pt_data/mupexi_cleavage_preds_iv.txt -o ./data/validation_data/pt_data/pepsickle_pt_summaries/mupexi_pt_summary.csv
