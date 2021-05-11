# Pepsickle paper repository
This repo contains the companion code for the manuscript: `Pepsickle rapidly and accurately predicts proteasomal cleavage sites for improved neoantigen identificat which can be found [here]().


## Trained model availability and use:
All fully trained models have been deployed as a separate software
package with instructions for installation and use. This can be found
at:
[https://github.com/pdxgx/pepsickle](https://github.com/pdxgx/pepsickle)

## How to run analysis in Linux/Unix:

1. Download the code in this repository:

   ```
   git clone https://github.com/pdxgx/pepsickle-paper.git
   cd ./pepsickle-paper
   ```

2. Setup and install necessary libraries (requires [python 3](https://www.python.org/downloads/) and [mysql](https://dev.mysql.com/doc/mysql-getting-started/en/#mysql-getting-started-installing) to be installed). For full list of python requirements see [reqirements.txt]():

   ```
   pip install -r requirements.txt
   pip install pepsickle
   ```

3. Enter the `data` directory and download the IEDB static data dump needed for analysis:

   NOTE: Paper analysis was performed on data pulled [INSERT DATE](). For identical reproduction subset all extracted data, including IEDB data, to entries on or before the specified date.
   ```
   cd ./data/raw/database_pulls
   wget http://www.iedb.org/downloader.php?file_name=doc/iedb_public.sql.gz
   ```
   The IEDB database ERD can also be found [here](http://www.iedb.org/downloader.php?file_name=doc/iedb_public_erd.pdf).
   
   The [AntiJen](http://www.ddg-pharmfac.net/antijen/AntiJen/aj_tcr.htm) database query feature is not working currently and has yet to be repaired, however the processed data from previously working queries can be found [here]().

4. Return to the main directory and edit the `MASTER.sh` script to include your mysql user name and
   password using nano or a text editor of your choice.

   `nano MASTER.sh`

5. Create output directories that are expected by the pipeline. The following directories are needed for proper pipeline output:
   ```
   mkdir XXXX [may just add this to .sh script w/ conditional
   ```

6. run the following command to iterate through data retrieval and
   processing steps:
   
   This script runs through the primary analysis and model training pipeline. Alternative models and options mentioned in the manuscript are also available but commented out for streamlining. Some steps are slow and annotated as such in comment lines.

   `bash run_analysis.sh`


## Information on data sources:

This project aggregates data from a variety of databases, including:
- [Immune Epitope Database (IEDB)](https://www.iedb.org/)
- [AntiJen database](http://www.ddg-pharmfac.net/antijen/AntiJen/antijenhomepage.htm)
- [SYFPEITHI database](http://www.syfpeithi.de/)

Data from peer reviewed literature was also aggregated. More details on paper specific data can be found under:
- `/docs/README - Digestion Sources.txt`
- `/docs/README - Winter et al data.txt`
- `/docs/README - Breast cancer epitope data.txt`

