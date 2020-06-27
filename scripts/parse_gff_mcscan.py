#!/usr/bin/env python3

import os, sys, re

gff = sys.argv[1]

with open(gff, 'r') as fh:
    for line in fh:
        if '#' in line:
            continue
        line = line.strip().split()
        if line[2] != 'gene':
            continue
        else:
            attrib = line[8].split(";")
            if 'genemodels' in gff:
                gene_id = [x for x in attrib if "ID=" in x][0]
                name = re.sub("ID=", "", gene_id) + "-mRNA-1"
            else:
                gene_id = [x for x in attrib if "Name=" in x][0]
                name = re.sub("Name=", "", gene_id)
            string_out = "{}\t{}\t{}\t{}".format(line[0], line[3], line[4], name)
            print(string_out)
