library(ape)
library(phytools)
library(tidyverse)

tree <- read.tree('genomes_tree/RAxML_bipartitionsBranchLabels.concatenated_alignments')
tr   <- root(tree, 5, r = TRUE)

metadata <- read_csv('tables/metadata.csv') %>%
  mutate(Agriculture = case_when(
    Agriculture == 'Outgroup1' ~ 'Hypomyces_Cladobotryum',
    Agriculture == 'Outgroup2' ~ 'Trichoderma',
    TRUE ~ Agriculture
  )) %>%
  mutate(acc = sub("\\..*$", "", acc))

id_map <- read_tsv('genomes_tree/id_map.txt', col_names = c('tree_id', 'acc')) %>%
  mutate(acc = sub("\\..*$", "", acc)) %>%
  mutate(acc = sub('GCF', 'GCA', acc)) %>%
  data.frame(., stringsAsFactors = FALSE) %>%
  inner_join(., metadata, by = 'acc') %>%
  mutate(complex_tree_id = paste0(Agriculture, ':', acc))

all_data <- id_map %>%
  select(acc, GenomeSize, complex_tree_id) %>%
    filter(!is.na(acc), 
         !is.na(GenomeSize))


rownames(id_map) <- id_map$tree_id
tr$tip.label     <- id_map[tr$tip.label, 'complex_tree_id']

states        <- all_data$GenomeSize
names(states) <- all_data$complex_tree_id

sub_tree <- keep.tip(tr, tip = names(states))
sub_tree <- rotate(sub_tree, 20)


fit.ace.true<-ace(states, sub_tree, type="continuous", method = "REML")

node_text <- paste0(signif(fit.ace.true$ace, 4), '\n(', signif(fit.ace.true$CI95[,1], 4), ' - ', signif(fit.ace.true$CI95[,2], 4), ')')

pdf(file = 'plots/ace_genome_size.pdf', width = 20, height = 10)
plotTree(sub_tree,mar=c(0.1,0.1,5.1,0.1),offset=0.5)
nodelabels(node_text, cex=0.7, bg = 'white')
dev.off()

### DO THE SAME WITH NUMBERS OF PROTEINS

tree <- read.tree('genomes_tree/RAxML_bipartitionsBranchLabels.concatenated_alignments')
tr   <- root(tree, 5, r = TRUE)

all_data <- id_map %>%
  filter(!is.na(acc), 
         !is.na(n_proteins)) %>%
  select(acc, n_proteins, complex_tree_id)

rownames(id_map) <- id_map$tree_id
tr$tip.label     <- id_map[tr$tip.label, 'complex_tree_id']


states        <- all_data$n_proteins
names(states) <- all_data$complex_tree_id

sub_tree <- keep.tip(tr, tip = names(states))

#plot(sub_tree)
#nodelabels()

sub_tree <- rotate(sub_tree, 37)
sub_tree <- rotate(sub_tree, 39)

fit.ace.true<-ace(states, sub_tree, type="continuous", method = "REML")

node_text <- paste0(signif(fit.ace.true$ace, 4), '\n(', signif(fit.ace.true$CI95[,1], 4), ' - ', signif(fit.ace.true$CI95[,2], 4), ')')

pdf(file = 'plots/ace_protein_number.pdf', width = 20, height = 10)
plotTree(sub_tree,mar=c(0.1,0.1,5.1,0.1),offset=0.5)
nodelabels(node_text, cex=0.7, bg = 'white')
dev.off()


