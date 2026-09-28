#!/usr/bin/env python3
import sys
input = open(sys.argv[1], 'r')
lines = input.readlines()[1:]

output_file = "00001.fasta"
for count, line in enumerate(lines, 1):
        # split into 25 peptides for each file
        if ((count - 1) % 25 == 0):
                num = "{0:0>5}".format(count)
                output_file = f"{num}.fasta"
        with open(f"{output_file}", 'a') as output:
                splitted = line.split() # Adjust for each input file 
                type=splitted[0]
                tissue=splitted[1]
                peptide = splitted[4]
                name= f"{peptide}_{tissue}_{type}_{count}"
                output.write("> ")
                output.write(name)
                output.write("\n")
                output.write(peptide)
                output.write("\n")
