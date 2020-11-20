library(ape)
library(tidyverse)
library(ggtree)
library(treeio)
library(ggplot2)
#tree_files <- list.files('gene_trees', full.names = TRUE, pattern = 'treefile')

tree_names <- scan('gene_trees/all_outgroups/gene_names.txt', what = 'character') 

trees <- read.tree(file = 'gene_trees/all_outgroups/in.tree', tree.names = tree_names)

uniq.trees <- unique.multiPhylo(trees)


tree <- read.newick('gene_trees/all_outgroups/in.tree', tree.names = tree_names)

ggdensitree(tree[1:100], alpha=.3, colour='steelblue') + 
  geom_tiplab(size=3) + xlim(0, 45)

ggtree(tree[1:100]) + facet_wrap(~.id, ncol=10)


btrees <- read.tree(system.file("extdata/RAxML", "RAxML_bootstrap.H3", package="treeio"))


ggdensitree(btrees, alpha=.3, colour='steelblue') + 
  geom_tiplab(size=3) + xlim(0, 45)
