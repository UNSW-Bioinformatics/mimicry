#!/bin/bash
indir=$1
file=$2
### These are specific to our HPC 
module load blast-plus/2.16.0
export BLASTDB=/data/bio/blastv5/

### Additional files
micr_file=ref/microbial.taxids

### Setting up inputs/outputs
base=$(basename ./input/${indir}/${file} .fasta)
qfile="./input/${indir}/${base}.fasta"
output="./output/${indir}/${base}.tsv"

echo $base
echo $qfile
echo $output

### Blast search
blastp -task blastp-short -query ${qfile} -db nr  -outfmt "6 qaccver qlen saccver slen sseq pident nident length mism
atch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames qcovs sskingdom sskingdoms stitle" -
taxidlist ${micr_file} -evalue 200000 -window_size 40 -threshold 11 -word_size 2 -gapopen 9 -gapextend 1 -matrix PAM3
0 -comp_based_stats 0  -out  ${output}

