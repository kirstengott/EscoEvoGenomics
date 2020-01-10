#!/usr/bin/env python3

import sys
import re
import pysam
import os

ortho = open(sys.argv[1], 'r')

## initialize a dictionary of dictionaries to store the orthologous terms
## I want to use this style of dictionary because I want to order one of my output tables by the column


ortho_dict = {}

tall_matrix = open('orthologues_full.csv', 'w')
wide_matrix = open('orthologues_presence_absence.csv', 'w')
prot_dir = 'proteins'

unique_genomes = []

def subset_fasta(file_in, file_out, sequences, genome):
    f_o = open(file_out, 'a')
    fasta_file = pysam.FastxFile(file_in)
    for entry in fasta_file:
        if entry.name in sequences:
            fa_name = entry.name + "__" + genome
            out = '>{0}\n{1}\n'.format(fa_name, entry.sequence)
            f_o.write(out)
        else:
            continue
    f_o.close()



## each line of the fastortho file is the ortholog followed by a list of genes and genomes
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
            ortho_dict[ortho_group][genome][0] += 1
            ortho_dict[ortho_group][genome][1].append(gene)
            
        else:
            ortho_dict[ortho_group][genome] = [1, [gene]]
tall_matrix.close()



genome_path_dict = {}

for i in unique_genomes:
    match = [x for x in os.listdir(prot_dir) if i in x][0]
    genome_path_dict[i] = match

    
title_string = "OrthoGroup,{}\n".format(",".join(unique_genomes))
n = 1


#genome_name = ['seq' + str(x) for x in range(0, len(unique_genomes))]
#genome_ids_out = dict(zip(unique_genomes, genome_name))

if not os.path.exists('orthologue_fastas'):
    os.mkdir('orthologue_fastas')

    
for ortho in ortho_dict:
    if n == 1:
        wide_matrix.write(title_string)
        n += 1
    orth_d = ortho_dict[ortho]
    ## return a list of the number of genes in all of the genomes
    #value_list = [str(orth_d.get(x, 0)) for x in unique_genomes]
    ## code to create a fasta file for each orthogroup
    ortho_fasta = os.path.join('orthologue_fastas/', ortho + ".fa")
    if os.path.exists(ortho_fasta):
        os.remove(ortho_fasta)
    value_list = []
    for x in unique_genomes:
        value_list.append(str(orth_d.get(x, 0)))
        if x in orth_d.keys():
            genes_list = orth_d[x][1]
            fasta = os.path.join('proteins', genome_path_dict[x])
            subset_fasta(file_in = fasta, file_out = ortho_fasta, sequences = genes_list, genome = x)
        else:
            pass
    value_o = ",".join(value_list)
    output_string = "{},{}\n".format(ortho, value_o)
    wide_matrix.write(output_string)

        
    
    
        
    


