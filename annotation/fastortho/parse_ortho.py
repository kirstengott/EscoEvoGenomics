#!/usr/bin/env python3

import sys
import re
import pysam
import os



## initialize a dictionary of dictionaries to store the orthologous terms
## I want to use this style of dictionary because I want to order one of my output tables by the column

def usage():
    usage = '''Usage: %s <FastOrthoOutput> <SequenceDirectory> <OrthoFastaOutdir> <SeqTypes>
               <FastOrthoOutput>: fastortho's output file
               <SequenceDirectory>: directory of sequences to write out orthogroups from
               <OrthoFastaOutdir>: where to output orthogroup fastas
               <SeqTypes>: the type of sequences used one of [cds, prot]\n''' % os.path.basename(sys.argv[0])
    sys.stderr.write(usage)
    sys.exit(1)


def main(argv):
    if len(argv) < 4:
        usage()
    elif "-h" in sys.argv[1]:
        usage()
    ortho = open(sys.argv[1], 'r')
    prot_dir = sys.argv[2]
    out_dir = sys.argv[3]
    seq_type = sys.argv[4]

    tall_matrix = open('orthologues_full.csv', 'w')
    wide_matrix = open('orthologues_presence_absence.csv', 'w')
    unique_genomes = []    
    ortho_dict = {} ## dictionary of dictionaries structured as {ortho_group : {genome : [genes]}}
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
            if '.all.maker.proteins' in genome:
                genome = re.sub('.all.maker.proteins', '', genome)
            elif '_protein'  in genome:
                genome = re.sub('_protein', '', genome)
            elif '.proteins' in genome:
                genome = re.sub('.proteins', '', genome)
            else:
                pass
            if genome not in unique_genomes:
                unique_genomes.append(genome)
            tall_matrix.write(",".join([ortho_group, genome, gene]) + "\n")
            if genome in ortho_dict[ortho_group]:
                ortho_dict[ortho_group][genome][0] += 1
                ortho_dict[ortho_group][genome][1].append(gene)
            else:
                ortho_dict[ortho_group][genome] = [1, [gene]]
    tall_matrix.close()

    ## parse in genomes
    genome_seqs_dict = {}
    for i in unique_genomes:
        match = [x for x in os.listdir(prot_dir) if i in x][0]
        fa_file = os.path.join(prot_dir, match)
        fasta_file = pysam.FastxFile(fa_file)
        genome_seqs_dict[i] = [fa_file, {}]
        seq1 = 0
        for entry in fasta_file:
            genome_seqs_dict[i][1][entry.name] = entry.sequence
            if seq1 == 0:
                seq_types = set()
                [seq_types.add(x) for x in entry.sequence]
                
        fasta_file.close()


    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    ## write out wide matrix and fasta files
    title_string = "OrthoGroup,{}\n".format(",".join(unique_genomes))
    n = 1
    for ortho in ortho_dict:
        if n == 1:
            wide_matrix.write(title_string)
            n += 1
        orth_d = ortho_dict[ortho]
        ## code to create a fasta file for each orthogroup
        ortho_fasta = os.path.join(out_dir, ortho + ".fa")
        if os.path.exists(ortho_fasta):
            os.remove(ortho_fasta)
        value_list = []
        for x in unique_genomes:
            if len(orth_d.keys()) > 3: ## only write out fasta's if there is more than 2 species represeted in the orthogroup
                value_list.append(str(orth_d.get(x, 0)))
                if x in orth_d.keys():
                    genes_list = orth_d[x][1]
                    fasta = genome_seqs_dict[x][0]
                    genes_dict = genome_seqs_dict[x][1]
                    f_o = open(ortho_fasta, 'a')
                    for gene in genes_list:
                        if seq_type == 'cds':
                            key = [x for x in genes_dict.keys() if gene in x][0]
                            out = '>{}_{}\n{}\n'.format(x, gene, genes_dict[key])
                        else:
                            out = '>{}_{}\n{}\n'.format(x, gene, genes_dict[gene])
                        f_o.write(out)
                    f_o.close()
                else:
                    pass
        value_o = ",".join(value_list)
        output_string = "{},{}\n".format(ortho, value_o)
        wide_matrix.write(output_string)
        
    

if __name__ == "__main__":
    main(sys.argv[1:])

