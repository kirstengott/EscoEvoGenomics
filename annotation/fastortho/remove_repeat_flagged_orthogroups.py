#!/usr/bin/env python3

import os
import sys

repeat_f_dir = '../repeat_flagged'
ortho_table = 'orthologues_full.csv'
ortho_directory = 'orthologue_fastas'



repeat_files = [os.path.join(repeat_f_dir, x) for x in os.listdir(repeat_f_dir)]
repeat_genes = []

for i in repeat_files:
    f = open(i, 'r')
    [repeat_genes.append(x.rstrip()) for x in f]
    f.close()




orth_table = open(ortho_table, 'r')


## here I am creating a dictionary or orthologues and the genes associated genes
ortho = {}

for line in orth_table:
    line = line.strip().split(",")
    orth = line[0]
    genome = line[1]
    gene = line[2]
    if orth not in ortho: 
        ortho[orth] = [gene]
    else:
        ortho[orth].append(gene)
orth_table.close()

ortho_flagged_out = open('repeat_flagged_orthogroups.txt', 'w')

num_orthogroups_flagged = 0
for orth in ortho:
    num_repeat_genes = len([x for x in ortho[orth] if x in repeat_genes])
    if num_repeat_genes > 0:
        total_genes = len(ortho[orth])
        if num_repeat_genes/total_genes >= 0.5:
            ortho_flagged_out.write(orth + "\n")
            orth_fa = orth + ".fa"
            fa_file = os.path.join(ortho_directory, orth_fa)
            os.remove(fa_file)
            num_orthogroups_flagged += 1

ortho_flagged_out.close()
print(num_orthogroups_flagged)
    
    
        
