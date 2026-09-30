### Scripts
This folder contains the different scripts used to generate the peptide candidates. As the underlying data takes time to generate, these scripts are provided to cover the broad steps, rather than to replicate the exact database. 

1. Extract peptides from CSV files into FASTA files. Note, this needs to be adjusted based on input data. -> csv2fasta.py
2. With the FASTA files, you will need to run each of the following scripts independently.
3. blast_search.sh -> performs the blast search using the FASTA file. Note, this assumes you have set up the nr database locally.
4. mhc_pred.py -> this runs mhcflurry on all the FASTA files in the provided directory. Note, this assumes you have installed MHCflurry and its associated tools, and installed its python library. 
5. netmhc_pred.sh -> this runs netMHCpan4.2 on each FASTA file. Note, this assumes you have this installed locally. Adjust script for your local setup.
6. extract_table.py -> this takes the blast outputs and overlaps them with the HMRGD data.
7. combine_hlas_pred.py -> this takes the HLA predictions and merges them with the match tables outputs. 

