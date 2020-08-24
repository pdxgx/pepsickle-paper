#!bin/bash
#### SETUP
# set working directory to base dir of project
cd /Users/weeder/PycharmProjects/pepsickle-paper
## set temp environmental vars for mysql use
export MYSQL_USER=[USER]
export MYSQL_PWD=[PASSWORD]

#### pull database data
## ANTIJEN
#### NOTE: issue with antijen database... see if temporary (query returns a table with no rows on web and via script)
#### python3 ./scripts/database_pulls/extract_AntiJen_data.py -o ./data/raw/database_pulls/
## SYFPEITHI
# python3 ./scripts/database_pulls/extract_SYF_data.py -a ./data/raw/raw_mammal_allele_series.csv -o ./data/raw/database_pulls/

#### extract static data
# IEDB
# wget -O ./data/raw/database_pulls/iedb_public.sql.gz https://www.iedb.org/downloader.php?file_name=doc/iedb_public.sql.gz
# gzip -d ./data/raw/database_pulls/iedb_public.sql.gz
# mysql --password=$MYSQL_PWD -u $MYSQL_USER < ./scripts/static_dataset_extractions/iedb_mysql_query.sql
## extract BC data
# python3 ./scripts/extract_breast_cancer_data.py -i ./data/raw/breast_cancer_data/ -o ./data/processed/
## extract digestion data
# python3 ./scripts/static_dataset_extractions/extract_digestion_data.py -i ./data/raw/digestion_map_files -o ./data/processed/

#### process extracted and static data
## database data
# AntiJen
# python3 ./scripts/dataset_processing/get_AntiJen_source_sequences.py -i ./data/raw/database_pulls/AntiJen_Tcell_epitopes.csv -o ./data/processed/
# SYFPEITHI
# python3 ./scripts/dataset_processing/get_SYFPEITHI_source_sequences.py -i ./data/raw/database_pulls/SYFPEITHI_epitopes.csv -o ./data/processed/

## static data
# IEDB
# python3 ./scripts/dataset_processing/iedb_mysql_processing.py -u $MYSQL_USER -p $MYSQL_PWD -o ./data/processed/
# BC data
# python3 ./scripts/dataset_processing/get_breast_cancer_source_sequences.py -i ./data/processed/breast_cancer_epitopes.csv -o ./data/processed/
# winter et al
# python3 ./scripts/dataset_processing/prep_winter_data.py -i ./data/raw/Winter_et_al_results.csv -o ./data/processed/
# levy et al
# python3 ./scripts/dataset_processing/prep_levy_data.py -i ./data/raw/Levy_et_al_degradome.csv -a ./data/raw/gene_protein_id_map.csv -o ./data/processed/levy_fragments_ambiguous.csv
# python3 ./scripts/merging_and_filtering/epitope_index_check.py -i ./data/processed/levy_fragments_ambiguous.csv -o ./data/processed/levy_fragments_ambiguous_mapped.csv
# python3 ./scripts/dataset_processing/levy_post_processing.py -i ./data/processed/levy_fragments_ambiguous_mapped.csv -o ./data/processed/levy_fragments_unique_mapped.csv

#### Merge datasets together
# python3 ./scripts/merging_and_filtering/merge_datasets.py -i ./data/processed -o ./data/merged
# NOTE run if only human data is desired
# python3 ./scripts/merging_and_filtering/merge_datasets.py -i ./data/processed -o ./data/merged --human-only

## Verify indices of epitopes in full file
# python3 ./scripts/merging_and_filtering/epitope_index_check.py -i ./data/merged/merged_data_all_mammal.csv -o ./data/merged/merged_data_all_mammal_clean.csv
# if only human data is desired
# python3 ./scripts/merging_and_filtering/epitope_index_check.py -i ./data/merged/merged_data_human_only.csv -o ./data/merged/merged_data_human_only_clean.csv

## split 20S and 26S data
# split for all mammal
# python3 ./scripts/merging_and_filtering/split_by_subunit_type.py -i ./data/merged/merged_data_all_mammal_clean.csv -o ./data/merged
# if only human
# python3 ./scripts/merging_and_filtering/split_by_subunit_type.py -i ./data/merged/merged_data_human_only_clean.csv -o ./data/merged --human-only

## generate the negative fragment examples based on annotated positives
### for all mammals
## for 20S
# python3 ./scripts/merging_and_filtering/negative_set_generation.py --full-negative-set -i ./data/merged/merged_20S_fragments.csv -o ./data/training_sets/all_mammal_20S_windows_13aa.pickle
## for 26S
# python3 ./scripts/merging_and_filtering/negative_set_generation.py --full-negative-set -i ./data/merged/merged_26S_fragments.csv -o ./data/training_sets/all_mammal_26S_windows_13aa.pickle
## for epitopes
# python3 ./scripts/merging_and_filtering/negative_set_generation.py --full-negative-set -i ./data/merged/merged_epitope_fragments.csv -o ./data/training_sets/all_mammal_epitope_windows_13aa.pickle

### for human only
## for 20S
# python3 ./scripts/merging_and_filtering/negative_set_generation.py --full-negative-set -i ./data/merged/merged_20S_fragments_human.csv -o ./data/training_sets/human_20S_windows_13aa.pickle
## for 26S
# python3 ./scripts/merging_and_filtering/negative_set_generation.py --full-negative-set -i ./data/merged/merged_26S_fragments_human.csv -o ./data/training_sets/human_26S_windows_13aa.pickle
## for epitopes
# python3 ./scripts/merging_and_filtering/negative_set_generation.py --full-negative-set -i ./data/merged/merged_epitope_fragments_human.csv -o ./data/training_sets/human_epitope_windows_13aa.pickle

## run model training
# 20S digestion
# python3 ./scripts/modeling/digestion_map_based_ensemble_net.py -i ./data/training_sets/all_mammal_20S_windows_13aa.pickle -o ./data/model_weights/
# 20S if using human only
# python3 ./scripts/modeling/digestion_map_based_ensemble_net.py --human-only -i ./data/training_sets/human_20S_windows_13aa.pickle -o ./data/model_weights/

# 26S


# epitope
# python3 ./scripts/modeling/epitope_based_ensemble_net.py -i ./data/training_sets/all_mammal_epitope_windows_13aa.pickle -o ./data/model_weights/
# epitope if using human only
# python3 ./scripts/modeling/epitope_based_ensemble_net.py --human-only -i ./data/training_sets/human_epitope_windows_13aa.pickle -o ./data/model_weights/

#### assess models on validation_prep data
## compile fragments for left out 20S digestion validation_prep data
# python3 ./scripts/static_dataset_extractions/extract_digestion_data.py -i ./data/validation_data/digestion_data/raw/ -o ./data/validation_data/digestion_data/
# python3 ./scripts/validation_prep/prep_20S_digestion_val_data.py -i ./data/validation_data/digestion_data/compiled_digestion_df.csv -o ./data/validation_data/digestion_data/
# python3 ./scripts/merging_and_filtering/epitope_index_check.py -i ./data/validation_data/digestion_data/20S_digestion_val_columns_remapped.csv -o ./data/validation_data/digestion_data/20S_digestion_val_data_indices_verified.csv

## compile epitope validation_prep data
# python3 ./scripts/validation_prep/prep_epitope_validation_data.py -i ./data/validation_data/epitope_data/41467_2016_BFncomms13404_MOESM1318_ESM.csv -o ./data/validation_data/epitope_data
# python3 ./scripts/merging_and_filtering/epitope_index_check.py -i ./data/validation_data/epitope_data/validation_epitopes_w_source.csv -o ./data/validation_data/epitope_data/validation_epitopes_verified.csv

## convert entries to cleavage windows
# python3 ./scripts/merging_and_filtering/negative_set_generation.py -i ./data/validation_data/digestion_data/20S_digestion_val_data_indices_verified.csv -o ./data/validation_data/validation_sets_pre-filter/20S_val_windows_13aa_paired.pickle
# python3 ./scripts/merging_and_filtering/negative_set_generation.py -i ./data/validation_data/epitope_data/validation_epitopes_verified.csv -o ./data/validation_data/validation_sets_pre-filter/epitope_val_windows_13aa_paired.pickle

## filter out entries seen in training data
# python3 ./scripts/validation_prep/filter_val_data.py -i ./data/validation_data/validation_sets_pre-filter -t ./data/training_sets -o ./data/validation_data/completed_validation_sets
