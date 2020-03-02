library(tidyverse)


metadata <- read_csv('tables/metadata.csv')


ggplot(metadata, aes(x = clade_groups, y = n_genes)) + 
  geom_boxplot() + 
  theme_minimal() +
  labs(x = 'Clade', 
       y = 'Number of Genes') +
  ggsave('plots/protein_coding_gene_N_boxplot.pdf')

ggplot(metadata, aes(x = clade_groups, y = assembly_length)) + 
  geom_boxplot() + 
  theme_minimal() +
  labs(x = 'Clade', 
       y = 'Assembly Length')  +
  ggsave('plots/assembly_length_boxplot.pdf')

ggplot(metadata, aes(x = clade_groups, y = GenomeSize)) + 
  geom_boxplot() + 
  theme_minimal() +
  labs(x = 'Clade', 
       y = 'Genome Size') +
  ggsave('plots/genome_size_boxplot.pdf')
