#!/usr/bin/env python

## this code was adapted from a networkx parser written by Marc Chevrette 
## this code adds by going back to the original bigscape input and outputing gff3 formatted files for the features in each component


import networkx as nx
import sys
import os
import re


import antismash2gff3






def main(argv):
    if len(argv) != 2:
        sys.stderr.write("Usage: %s <in.network> <bigscape_input_dir>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)

    edge_thresh = 0.3
    filename = sys.argv[1]
    bigscape_input = sys.argv[2]


    all_nodes = {}
    g = nx.Graph()

    outdir = 'gff3'
    try:
        os.mkdir(outdir)
    except OSError:
        print ("Creation of the directory %s failed, it may already exist." % outdir)
    else:
        print ("Successfully created the directory %s " % outdir)

    parsed_network_dir = 'parsed_networks'    
    try:
        os.mkdir(parsed_network_dir)
    except OSError:
        print ("Creation of the directory %s failed, it may already exist." % parsed_network_dir)
    else:
        print ("Successfully created the directory %s " % parsed_network_dir)


        
    bigscape_input_files = os.listdir(bigscape_input)
    bgsc_paths = [os.path.join(bigscape_input, x) for x in bigscape_input_files]


    file_root = os.path.basename(filename)

    with open(filename, 'r') as f:
            next(f)
            for ln in f:
                    a = ln.split('\t')
                    if a[0] not in all_nodes:
                            g.add_node(a[0])
                            all_nodes[a[0]] = 1
                    if a[1] not in all_nodes:
                            g.add_node(a[1])
                            all_nodes[a[1]] = 1
                    if a[0] != a[1]:
                            g.add_edge(a[0], a[1], weight = a[7])

    sg = nx.connected_component_subgraphs(g)

    c = 1
    bgc = ''
    parsed_network_outfile = parsed_network_dir + "/" + file_root
    parsed_n = open(parsed_network_outfile, 'a')
    for comp in sg:
            component = 'component_' + str(c)
            outfile = 'gff3/' + file_root + "_" + component + ".gff3"
            if len(comp) > 1: ## I don't care about singletons for this analysis
                    out_f = open(outfile, 'a')
                    all_gff3_lines = []
                    for node in comp:
                            parsed_n.write('\t'.join(['component_' + str(c), node]) + "\n")
                            if 'BGC' in node:
                                    bgc = "BGCid:" + node
                            else:
                                    matching = [s for s in bgsc_paths if node in s]
                                    if matching:
                                            gff3_lines = antismash2gff3.parse_antismash4(in_file = matching[0], cluster_id = node)
                                            all_gff3_lines.extend(gff3_lines)
                    if len(bgc) > 0:
                        all_gff3_lines = ["{0},{1}".format(x, bgc) for x in all_gff3_lines]
                    all_lines = "\n".join(all_gff3_lines)
                    bgc = ''
                    out_f.write(all_lines)
                    out_f.close()
            else:
                    pass


            c += 1
    parsed_n.close()

if __name__ == '__main__':
    main(sys.argv[1:])
