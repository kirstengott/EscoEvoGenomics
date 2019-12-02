#!/usr/bin/env python3

import sys
import re
import os

in_file = open(sys.argv[1], 'r')

f_out = os.path.basename(os.path.splitext(sys.argv[1])[0]) + "_repeat_flagged.txt"

flagged = open(f_out, 'w')

gene_dict = {}
for line in in_file:
    liner = line.strip().split("\t")
    ## get the gene name out
    gene_r = [x for x in liner[8].split(";") if "Parent" in x][0]
    gene = re.sub("Parent=", "", gene_r)
    repeat_overlap = int(liner[12])
    exon_length = abs(int(liner[4]) - int(liner[3]))
    if gene not in gene_dict:
        gene_dict[gene] = {'repeat' : repeat_overlap,
                           'length' : exon_length}
    else:
        gene_dict[gene]['repeat'] += repeat_overlap
        gene_dict[gene]['length'] += exon_length


        
for gene in gene_dict:
    perc_repeat = gene_dict[gene]['repeat'] / gene_dict[gene]['length']
    if perc_repeat >= 0.2:
        flagged.write(gene + "\n")
    else:
        continue
    
