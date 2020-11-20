#!/usr/bin/env python3

import sys, re
from collections import defaultdict
import glob
import pandas as pd

def parse_interproscan(interpro, genome_id):
    annots = defaultdict(dict)
    with open(interpro, 'r') as f:
        for line1 in f:
            line = line1.strip().split()
            gene = line[0]
            tool = line[3]
            annot = line[4]
            if tool == 'TMHMM':
                if not gene in annots[genome_id]:
                    annots[genome_id][gene] = {tool : 1}
                else:
                    if not tool in annots[genome_id][gene]:
                        annots[genome_id][gene][tool] = 1
                    else:
                        annots[genome_id][gene][tool] += 1
            else:
                if not gene in annots[genome_id]:
                    annots[genome_id][gene] = {tool : [annot]}
                else:
                    if not tool in annots[genome_id][gene]:
                        annots[genome_id][gene][tool] = [annot]
                    else:
                        annots[genome_id][gene][tool].append(annot)
                    
            ## parse out go terms
            if 'GO:' in line1:
                tool = 'GO'
                annot = [x for x in line if 'GO:' in x][0].split("|")
                if not tool in annots[genome_id][gene]:
                    annots[genome_id][gene][tool] = annot
                else:
                    for x in annot:
                        if not x in annots[genome_id][gene][tool]: ## only add in unique GO terms
                            annots[genome_id][gene][tool].append(x)
                        else:
                            continue
    return(annots)

def parse_blast(annots, file_in, genome_id, tool):
    with open(file_in, 'r') as f:
        for line in f:
            line = line.strip().split()
            gene = line[0]
            target = line[1]
            if not gene in annots[genome_id]:
                annots[genome_id][gene] = {tool : target}
            else:
                if not tool in annots[genome_id][gene]:
                    annots[genome_id][gene][tool] = target
                else:
                    annots[genome_id][gene][tool] += "," + target
        return(annots)

def parse_blast_merops(annots, file_in, genome_id, tool):
    with open(file_in, 'r') as f:
        for line in f:
            line = line.strip().split()
            gene = line[0]
            target = line[1]
            if not gene in annots[genome_id]:
                annots[genome_id][gene] = {tool : target}
            else:
                if not tool in annots[genome_id][gene]:
                    annots[genome_id][gene][tool] = target
                else:
                    annots[genome_id][gene][tool] += "," + target
        return(annots)

def parse_card(annots, card, genome_id):
    header = 0
    with open(card, 'r') as f:
        for line in f:
            if header == 0:
                header += 1
                continue
            else:
                line = line.strip().split("\t")
                count = 0
                gene = line[0].split()[0]
                drug = line[14]
                res = line[16] ## keep the gene family ARO 
                if not gene in annots[genome_id]:
                    annots[genome_id][gene] = {'CARD' : res}
                else:
                    if not 'CARD' in annots[genome_id][gene]:
                        annots[genome_id][gene]['CARD'] = res
                    else:
                        annots[genome_id][gene]['CARD'].append(res)
        return(annots)
    
def parse_targetp(annots, targetp, genome_id):
    with open(targetp, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            line = line.strip().split()
            gene = line[0]
            loc = line[1]
            if loc == 'noTP':
                continue
            else:
                if not gene in annots[genome_id]:
                    annots[genome_id][gene] = {'targetp' : loc}
                else:
                    if not 'targetp' in annots[genome_id][gene]:
                        annots[genome_id][gene]['targetp'] = loc
                    else:
                        annots[genome_id][gene]['targetp'].append(loc)
        return(annots)


def parse_wolfpsort(annots, wolfpsort, genome_id):
    with open(wolfpsort, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split(",")
            gene = line[0].split()[0]
            loc1 = line[0].split()[1]
            loc1_score = line[0].split()[2]
            loc_out = [loc1]
            if len(line) > 1:
                loc2 = line[1].split()[0]
                loc2_score = line[1].split()[1]
                if float(loc1_score)/float(loc2_score) < 2:
                    loc_out.append(loc2)
            if not gene in annots[genome_id]:
                annots[genome_id][gene] = {'wolfpsort' : loc_out}
            else:
                if not 'wolfpsort' in annots[genome_id][gene]:
                    annots[genome_id][gene]['wolfpsort'] = loc_out
                else:
                    annots[genome_id][gene]['wolfpsort'].append(loc_out)
        return(annots)


def parse_cazy(cazy, annots, genome_id):
    header = 0
    with open(cazy, 'r') as f:
        for line in f:
            if header == 0:
                header +=1
                continue
            else:
                line = line.strip().split()
                gene = line[0]
                annot = line[1].split("(")[0]
                if not gene in annots[genome_id]:
                    annots[genome_id][gene] = {'cazy' : [annot]}
                else:
                    if not 'cazy' in annots[genome_id][gene]:
                        annots[genome_id][gene]['cazy'] = [annot]
                    else:
                        annots[genome_id][gene]['cazy'].append(annot)
    return(annots)


def parse_led(led, annots, genome_id):
    with open(led, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            else:
                line = line.strip().split()
                gene = line[0]
                annot = re.sub("\..*$", "", line[1])
                if not gene in annots[genome_id]:
                    annots[genome_id][gene] = {'led' : [annot]}
                else:
                    if not 'led' in annots[genome_id][gene]:
                        annots[genome_id][gene]['led'] = [annot]
                    else:
                        if annot in annots[genome_id][gene]['led']:
                            continue
                        else:
                            annots[genome_id][gene]['led'].append(annot)
    return(annots)

def parse_busco(busco, annots, genome_id):
    with open(busco, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            else:
                if 'Missing' in line:
                    continue
                line = line.strip().split()
                gene = line[2]
                annot = " ".join(line[6:])
                if not gene in annots[genome_id]:
                    annots[genome_id][gene] = {'busco' : [annot]}
                else:
                    if not 'busco' in annots[genome_id][gene]:
                        annots[genome_id][gene]['busco'] = [annot]
                    else:
                        if annot in annots[genome_id][gene]['busco']:
                            continue
                        else:
                            annots[genome_id][gene]['busco'].append(annot)
    return(annots)



genomes = ['SPDT00000000',
           'GCA_000167675',
           'GCA_000170995',
           'GCA_000171015',
           'GCA_000513815',
           'GCA_000988865',
           'GCA_001050175',
           'GCA_001481775',
           'GCA_002022785',
           'GCA_002838845',
           'GCA_002894145',
           'GCA_002894205',
           'GCA_003012105',
           'GCA_003025095',
           'GCA_003025105',
           'GCA_003025115',
           'GCA_003025155',
           'GCA_011066345',
           'GCA_004303015',
           'ICBG1054',
           'ICBG1065',
           'ICBG1075',
           'ICBG1096',
           'ICBG2046',
           'ICBG2047',
           'ICBG2048',
           'ICBG2049',
           'ICBG710',
           'ICBG712',
           'ICBG721',
           'ICBG726',
           'ICBG730',
           'ICBG731',
           'ICBG733',
           'ICBG736',
           'ICBG742',
           'ICBG751',
           'LGSR00000000',
           'NIGB00000000',
           'NIGC00000000',
           'NIGD00000000',
           'NQYQ00000000',
           'NQYR00000000',
           'NQYS00000000',
           'GCA_011799845.1']




#orthos = 'Orthogroups.tsv'



## have to read in with pandas because the whitespace is weird
#df = pd.read_csv(orthos, sep="\t")
#df.dropna(thresh = 1)
#orth    = df.set_index('Orthogroup').T.to_dict() ## make it a dictionary



# orth = defaultdict(dict)

# with open(orthos) as f:
#     for line in f:
#         all_l = line.strip().split(",")
#         orthogroup = all_l[0]
#         genome_path = all_l[1]
#         genome_id = [x for x in genomes if x in genome_path][0]
#         gene = all_l[2]
#         if not genome_id in orth[orthogroup]:
#             orth[orthogroup][genome_id] = [gene]
#         else:
#             orth[orthogroup][genome_id].append(gene)


all_annots_out = open('all_annotations.txt', 'w')

for genome_id in genomes:

    files = glob.glob("*/*" + genome_id + "*")
    ## get file names
    cazy = [x for x in files if 'cazy' in x][0]
    interpro  = [x for x in files if 'interpro' in x][0]
    led       = [x for x in files if 'LED' in x][0]
    merops    = [x for x in files if 'merops' in x][0]
    card      = [x for x in files if 'card' in x][0]
    targetp   = [x for x in files if 'targetp' in x][0]
    wolfpsort = [x for x in files if 'wolfpsort' in x][0]
    virulence = [x for x in files if 'virulence' in x][0]
    busco     = [x for x in files if 'ascomybocta_odb10_busco' in x][0]
    
    ## get annotations
    annots = parse_interproscan(interpro = interpro, genome_id = genome_id)
    annots = parse_cazy(cazy = cazy, annots = annots, genome_id = genome_id)
    annots = parse_led(annots = annots, led = led, genome_id = genome_id)
    annots = parse_blast(annots = annots, file_in = merops, genome_id = genome_id, tool = 'MEROPS')
    annots = parse_card(annots = annots, card = card, genome_id = genome_id)
    annots = parse_targetp(annots = annots, targetp = targetp, genome_id = genome_id)
    annots = parse_wolfpsort(annots = annots, wolfpsort = wolfpsort, genome_id = genome_id)
    annots = parse_blast(annots = annots, file_in = virulence, genome_id = genome_id, tool = 'virulence')
    annots = parse_busco(annots = annots, busco = busco, genome_id = genome_id)

    for gene in annots[genome_id]:
        sub = annots[genome_id][gene]
        for tool in sub:
            if isinstance(sub[tool], list):
                if len(sub[tool]) > 1:
                    ann = ",".join(sub[tool])
                else:
                    ann = sub[tool][0]
            else:
                ann = sub[tool]
            all_annots_out.write("{}\t{}\t{}\t{}\n".format(genome_id, gene, tool, ann))

all_annots_out.close()



out_f = open('all_annotations.txt', 'a')

header = 0
with open('antismash_genes.csv', 'r') as fh:
    for line in fh:
        if header == 0:
            header += 1
            continue
        else:
            line = line.strip().split(",")
            genome = line[2]
            gene   = line[1]
            cluster_type = line[0]
            out_str = "{}\t{}\t{}\t{}\n".format(genome, gene, 'antismash', cluster_type)
            out_f.write(out_str)
out_f.close()
