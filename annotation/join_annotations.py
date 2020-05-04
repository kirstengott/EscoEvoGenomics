#!/usr/bin/env python3

import sys
from collections import defaultdict
import glob

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
            if not gene in annots[genome_id]:
                annots[genome_id][gene] = {tool : 1}
            else:
                if not tool in annots[genome_id][gene]:
                    annots[genome_id][gene][tool] =  1
                else:
                    annots[genome_id][gene][tool] += 1
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
                res = line[15]
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
            line = line.strip().split()
            gene = line[0]
            loc = line[5]
            if loc == '_':
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

    
genomes = ['SPDT00000000',
           'GCA_000167675',
          'GCA_000170995',
          'GCA_000171015',
          'GCA_001050175',
          'GCA_001481775',
          'GCA_002022785',
          'GCA_003012105',
          'GCA_003025095',
          'GCA_003025105',
          'GCA_003025115',
          'GCA_003025155',
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
           'NQYS00000000']



orthos = 'fastortho/orthologues_full.csv'

orth = defaultdict(dict)

with open(orthos) as f:
    for line in f:
        all_l = line.strip().split(",")
        orthogroup = all_l[0]
        genome_path = all_l[1]
        genome_id = [x for x in genomes if x in genome_path][0]
        gene = all_l[2]
        if not genome_id in orth[orthogroup]:
            orth[orthogroup][genome_id] = [gene]
        else:
            orth[orthogroup][genome_id].append(gene)


all_annots_out = open('all_annotations.txt', 'w')

for genome_id in genomes:

    

    files = glob.glob("*/*" + genome_id + "*")
    print(genome_id)

    
    ## get file names
    cazy = [x for x in files if 'cazy' in x][0]
    interpro  = [x for x in files if 'interpro' in x][0]
    led       = [x for x in files if 'LED' in x][0]
    merops    = [x for x in files if 'merops' in x][0]
    card      = [x for x in files if 'card' in x][0]
    targetp   = [x for x in files if 'targetp' in x][0]
    wolfpsort = [x for x in files if 'wolfpsort' in x][0]
    virulence = [x for x in files if 'virulence' in x][0]
    
    ## get annotations
    annots = parse_interproscan(interpro = interpro, genome_id = genome_id)
    annots = parse_cazy(cazy = cazy, annots = annots, genome_id = genome_id)
    annots = parse_blast(annots = annots, file_in = led, genome_id = genome_id, tool = 'LED')
    annots = parse_blast(annots = annots, file_in = merops, genome_id = genome_id, tool = 'MEROPS')
    annots = parse_card(annots = annots, card = card, genome_id = genome_id)
    annots = parse_targetp(annots = annots, targetp = targetp, genome_id = genome_id)
    annots = parse_wolfpsort(annots = annots, wolfpsort = wolfpsort, genome_id = genome_id)
    annots = parse_blast(annots = annots, file_in = virulence, genome_id = genome_id, tool = 'virulence')

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


orth_annot = 'orthogroups_annotation.txt'

o_annot_out = open(orth_annot, 'w')

for o in orth:
    genomes = orth[o].keys()
    annot = {}
    for genome in genomes:
        genes = orth[o][genome]
        for gene in genes:
            try:
                for an in annots[genome][gene]:
                    if not an in annot:
                        if not isinstance(annots[genome][gene][an], list):
                            annot[an] = [annots[genome][gene][an]]
                        else:
                            if not an in annot:
                                annot[an] = annots[genome][gene][an]
                    else:
                        if not isinstance(annots[genome][gene][an], list):
                            if not annots[genome][gene][an] in annot[an]:
                                annot[an].append(annots[genome][gene][an])
                            else:
                                continue
                        else:
                            for sub in annots[genome][gene][an]:
                                if not sub in annot[an]:
                                    annot[an].append(sub)
                                else:
                                    continue
            except:
                continue
    for x in annot:
        if len(annot[x]) > 1:
            join_annot = [str(x) for x in annot[x]]
            ann = ",".join(join_annot)
        else:
            ann = str(annot[x][0])
        o_annot_out.write(("{ortho} {genomes} {tool} {annot}\n".format(ortho = o, genomes = ",".join(genomes), tool = x, annot = ann)))
    
o_annot_out.close()

