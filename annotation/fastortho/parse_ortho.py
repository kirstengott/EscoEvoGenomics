#!/usr/bin/env python3

import sys
import re

ortho = open(sys.argv[1], 'r')

## initialize a dictionary of dictionaries to store the orthologous terms
## I want to use this style of dictionary because I want to order one of my output tables by the column


ortho_dict = {}

tall_matrix = open('orthologues_full.csv', 'w')
wide_matrix = open('orthologues_presence_absence.csv', 'w')

unique_genomes = []

for line in ortho:
    line = line.strip().split(':')
    ortho_group = line[0].split(' ')[0]
    genes_list = line[1].strip().split(" ")
    ortho_dict[ortho_group] = {}
    for g in genes_list:
        gs = g.split("(")
        gene = gs[0]
        genome = re.sub("\)", "", gs[1])
        if genome not in unique_genomes:
            unique_genomes.append(genome)
        tall_matrix.write(",".join([ortho_group, genome, gene]) + "\n")
        if genome in ortho_dict[ortho_group]:
            ortho_dict[ortho_group][genome] += 1
        else:
            ortho_dict[ortho_group][genome] = 1
tall_matrix.close()




title_string = "OrthoGroup,{}\n".format(",".join(unique_genomes))
n = 1
for ortho in ortho_dict:
    if n == 1:
        wide_matrix.write(title_string)
        n += 1
    orth_d = ortho_dict[ortho]
    value_list = [str(orth_d.get(x, 0)) for x in unique_genomes]
    value_o = ",".join(value_list)
    output_string = "{},{}\n".format(ortho, value_o)
    wide_matrix.write(output_string)
    
        
    

