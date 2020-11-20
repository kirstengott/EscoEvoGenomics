#!/usr/bin/env python3

import sys, os, re
from Bio import SeqIO, Seq

out_f = open(sys.argv[3], 'w')
record_dict = SeqIO.index(sys.argv[2], "fasta")

gff_cds = {}

with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        line = line.strip().split()
        chrom = line[0]
        name = line[1]
        start = int(line[2])
        stop = int(line[3])
        if start < stop:
            strand = '+'
        else:
            strand = '-'
        gff_cds[name] = [chrom, start, stop, strand]
            
## goes through the input fasta file and pulls out CDS sequences

for mrna in gff_cds:
    cds = gff_cds[mrna]
    start  = cds[1]
    stop   = cds[2] ## grab the last stop and strand
    strand = cds[3]
    
    seq = record_dict[cds[0]][start-1:stop].seq
    if strand == "-":
        seq = seq.reverse_complement()
    while len(seq) % 3 > 0: ## adding N's to partial CDS
        seq += 'N'
    protein = seq.translate()
    outstring = ">{m_id} # {start} # {stop} # {strand}\n{seq}\n".format(m_id = mrna, start = start, stop = stop, strand = strand, seq = protein)
    out_f.write(outstring)
out_f.close()
record_dict.close()
