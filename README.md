# Pepsickle Paper

Modeling and prediction of proteasomal cleavage sites. Find the full
paper at: `[ref]`

## Trained model availability and use:

All fully trained models have been deployed as a separate software
package with instructions for installation and use. They can be found
at: [https://github.com/pdxgx/pepsickle](https://github.com/pdxgx/pepsickle)

## How to run analysis in Linux/Unix:

1. Download the code in this repository:

   `git clone https://github.com/pdxgx/pepsickle-paper.git`

2. Setup and install necessary libraries (requires python 3 to be
   installed):

   `pip install .`

3. Enter the base directory Download static data files needed:

   `wget [path]`

4. Edit the `run_analysis` script to include your mysql user name and
   password.

5. run the following command:

   `bash run_analysis.sh`


## Information on data sources:
This project aggregates data from a variety of databases, including:
- [Immune Epitope Database (IEDB)](https://www.iedb.org/)
- [AntiJen database](http://www.ddg-pharmfac.net/antijen/AntiJen/antijenhomepage.htm)
- [SYFPEITHI database](http://www.syfpeithi.de/)

Data from peer reviewed literature was also aggregated. More details on paper specific details and references can be found under:
- `/docs/README - Digestion Sources.txt`
- `/docs/README - Winter et al data.txt`
- `/docs/README - Breast cancer epitope data.txt`
