import os
import numpy as np
from Bio import SeqIO
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re
from mhcflurry import Class1AffinityPredictor

alleles = ["HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:03", "HLA-A*02:06", "HLA-A*03:01",
        "HLA-A*11:01", "HLA-A*23:01", "HLA-A*24:02", "HLA-A*26:01", "HLA-A*30:01",
        "HLA-A*30:02", "HLA-A*31:01", "HLA-A*32:01", "HLA-A*33:01", "HLA-A*68:01",
        "HLA-A*68:02", "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*35:01",
        "HLA-B*40:01", "HLA-B*44:02", "HLA-B*44:03", "HLA-B*51:01", "HLA-B*53:01",
        "HLA-B*57:01", "HLA-B*58:01"]

mhcflurry_threshold=100.0

def read_in_peptides(fasta_file):
    peptides = {} 
    for seq_record in SeqIO.parse(open(fasta_file, mode='r'), 'fasta'):
        id = seq_record.id
        seqs = seq_record.seq
        pep_test = str(seqs)
        pep_test = pep_test.replace(r'*', "")
        peptides[id] = pep_test
    peptides = list(peptides.values() ) 
    return peptides

def run_mhcflurry(peptides, alleles, fileout):
    handle = open(fileout, "w")
    predictor = Class1AffinityPredictor.load()
    preds = {p: {} for p in peptides}
    for allele in alleles:
        try:
            df = predictor.predict_to_dataframe(
                peptides=peptides,
                allele=allele,
                include_percentile_ranks=True,
                include_confidence_intervals=False,
                throw=True
            )
            for _, row in df.iterrows():
                pep = row['peptide']
                rank = row['prediction_percentile']
                print(row['peptide'], row['allele'], row['prediction'], row['prediction_percentile'], file=handle)
                preds[pep][f"MHCflurry_{allele}"] = rank
        except Exception as e:
            print(f"Failed MHCflurry for {allele}: {str(e)}")
    handle.close() 
    return preds

## Use folder with fasta files 
fasta_dir = "./" 
fasta_bases = [f.rpartition('.')[0] for f in os.listdir(fasta_dir) if re.match(r'[0-9]+.*\.fasta', f)]
for fasta_base in fasta_bases[0:311]: 
    print(fasta_base)
    fasta_file = fasta_dir + fasta_base + ".fasta"
    peptides =  read_in_peptides(fasta_file)
    peptides_short = list(filter(lambda x: len(x) < 16, peptides)) 
    peptides_long = list(filter(lambda x: len(x) > 15, peptides)) 
    peptides_short = [x for x in peptides_short if not re.search(r'[BXZOU]+', x)]
    peptides_short = list(set(peptides_short))
    fasta_out= fasta_dir + "mhcflurry_" + fasta_base + ".out"
    mf = run_mhcflurry(peptides_short, alleles, fasta_out)


