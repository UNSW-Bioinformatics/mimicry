#!/bin/bash
indir=$1
file=$2 

### Specific to our HPC
export NETMHCpan=bin/netmhc/netMHCpan-4.2/Linux_x86_64/

### Set up inputs and outputs
base=$(basename ./input/${indir}/${file} .fasta)
qfile="./input/${indir}/${base}.fasta"
pfile="./input/${indir}/${base}_pep"
output="./hla_pred_out/${indir}/netMHC_${base}.out"

grep '>' ${qfile} -v > $pfile

$NETMHCpan/bin/netMHCpan-4.2 \
 -a HLA-A01:01,HLA-A02:01,HLA-A02:03,HLA-A02:06,HLA-A03:01,HLA-A11:01,HLA-A23:01,HLA-A24:02,HLA-A26:01,HLA-A30:01,HLA
-A30:02,HLA-A31:01,HLA-A32:01,HLA-A33:01,HLA-A68:01,HLA-A68:02,HLA-B07:02,HLA-B08:01,HLA-B15:01,HLA-B35:01,HLA-B40:01
,HLA-B44:02,HLA-B44:03,HLA-B51:01,HLA-B53:01,HLA-B57:01,HLA-B58:01 \
 -f ${pfile} -p -BA > ${output}
