library(tidyverse)
library(ggtree)
library(ape)

load('rdata/orthologues_long.rdata')
load('rdata/orthologues_long_annot.rdata')


orthog_long %>% filter(Orthogroup %in% genes_include) %>% 
  .$genes %>% 
  unique() %>% 
  write(., file = 'pangenome_genes.txt')

genes_include = scan('anvio_genes_include.txt', what = 'character')

gene_network <- orthog_long %>%
  group_by(Orthogroup, acc) %>%
  summarize(num_genes = length(genes)) %>%
  ungroup() %>% filter(Orthogroup %in% genes_include)

## for anvio

meta <- read_csv('tables/metadata.csv') 

tree <- read.tree(file = 'astral_tree/out_tree-scored-t32.tre')

tree$node.label <- tree$node.label %>% sub("\\[.*$", "", .)
tree_sub = extract.clade(tree, node = 'N68')



## load object 'genomes_keep'
source("scripts/working_genomes.R")

#cat(paste0('"', tree_sub$tip.label[which(tree_sub$tip.label %in% genomes_keep)], '"'), sep = ",\n")



############# ANVIO TREE LAYER #########################
tree_out <- drop.tip(tree_sub, tip = tree_sub$tip.label[!tree_sub$tip.label %in% genomes_keep])

write.tree(tree_out, 'anvio/layers-tree.nwk')

########################################################




################### ANVIO DATA LAYER ########################
### Orthologue PRESENCE OR ABSENCE DATA
t <- gene_network %>% 
  mutate(acc = ifelse(acc == "SPDT00000000", yes = "GCA_008477525", no = acc)) %>%
  mutate(acc = ifelse(grepl('GCA', acc), yes = paste0(acc, '.1'), no = acc)) %>% 
  mutate(presence = 1) %>% 
  dplyr::select(-num_genes) %>%
  filter(acc %in% genomes_keep) %>%
  spread(acc, presence) 
t[is.na(t)] <- 0

t %>% write_tsv('anvio/data.txt')

##################################################################




################ ANVIO ADDITIONAL ANNOTATIONS ##############################
### per orthologue

in_upset_plot <- scan('rdata/anvio_ortho.txt', what = 'character')

annots_keep <- c('Pfam', 'busco', 'virulence', 'MEROPS', 'CARD', 'cazy', 'led', 'antismash')

orthog_long_annot %>% 
  select(-genome) %>% 
  mutate(acc = ifelse(acc == "SPDT00000000", yes = "GCA_008477525", no = acc)) %>%
  mutate(acc = ifelse(grepl('GCA', acc), yes = paste0(acc, '.1'), no = acc)) %>% 
  filter(acc %in% genomes_keep) %>%
  distinct() %>%
  group_by(Orthogroup) %>%
  mutate(n_genomes = length(unique(acc))) %>%
  ungroup() %>%
  mutate(tools_keep = ifelse(tool %in% annots_keep, yes = tool, no = 'other')) %>%
  select(Orthogroup, tools_keep, n_genomes) %>%
  distinct() %>%
  mutate(Presence = 'True') %>%
  spread(tools_keep, Presence) %>%
  select(-other) %>%
  mutate_at(c('antismash', 'busco', 'CARD', 'cazy', 'led', 'MEROPS', 'virulence'), 
            function(x){ifelse(is.na(x), yes = 'False', no = x)}) %>% 
  filter(!is.na(n_genomes)) %>% 
  select(-n_genomes) %>%
#  mutate(in_upset_plot = ifelse(Orthogroup %in% in_upset_plot, yes = 'True', no = 'False')) %>%
  write_tsv('anvio/additional-items-data.txt')

############################################################################################


######################## ANVIO ANI DATA ######################################## (not used)


ani <- read_tsv('tables/ani.txt', col_names = c('Query', 'Reference', 'ANI', 't', 't1')) %>% 
  select(-t, -t1) %>%
  mutate_at(c('Query', 'Reference'), basename) %>%
  mutate_at(c('Query', 'Reference'), function(x){sub(".cds_from_genomic.fa", "", x)}) %>%
  mutate_at(c('Query', 'Reference'), function(x){sub("\\..*$", "", x)}) %>%
#  select(-QueryAligned) %>% 
  rename('acc' = Query) %>% 
  mutate(acc = ifelse(acc == "SPDT00000000", yes = "GCA_008477525", no = acc)) %>%
  mutate(Reference = ifelse(Reference == "SPDT00000000", yes = "GCA_008477525", no = Reference)) %>%
  mutate(acc = ifelse(grepl('GCA', acc), yes = paste0(acc, '.1'), no = acc)) %>%
  mutate(Reference = ifelse(grepl('GCA', Reference), yes = paste0(Reference, '.1'), no = Reference)) %>%
  filter(acc %in% genomes_keep) %>%
  spread(Reference, ANI) 

#######################################################################################

########## ANVIO ADDITIONAL LAYERS ######################################################
## per genome

meta  %>% 
  dplyr::select(acc, genus) %>%
  #  left_join(., ani, by = 'acc') %>%
  filter(acc %in% genomes_keep) %>% 
  write_tsv('anvio/additional-layers-data.txt')













############# DATA FOR ANVIO SIMPLE

ani_groups








