#!/bin/env python
​
import networkx as nx
​
edge_thresh = 0.3
filename = 'your_network.tsv'
​
all_nodes = {}
g = nx.Graph()
​
with open(filename, 'r') as f:
	next(f) ## header
	for ln in f.read().splitlines():
		a = ln.split("\t")
		if a[0] not in all_nodes:
			g.add_node(a[0])
			all_nodes[a[0]] = 1
		if a[1] not in all_nodes:
			g.add_node(a[1])
			all_nodes[a[1]] = 1
		if a[0] != a[1]:
			g.add_edge(a[0], a[1], weight = a[7]) ## a[7] is rawdistance
​
sg = nx.connected_component_subgraphs(g)
c = 1 ## iterator for component number
for comp in sg:
	for node in comp:
		print("\t".join(['component_'+str(c), node]))
	c += 1
