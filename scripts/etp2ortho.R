library(tidyverse)

cds2protein <- read_tsv('fasta/cds2protein.map', col_names = c('gene', 'protein'))

orthos <- read_tsv('antismash_cluster_analysis/ETP/cds_regions_focal_genome.orthogroups') %>%
  gather(genome, protein, -Orthogroup) %>%
  filter(!is.na(protein)) %>%
  mutate(protein = strsplit(protein, split = ", ")) %>%
  unnest(protein) %>%
  mutate(protein = gsub(" ", "", protein)) %>%
  left_join(cds2protein, by = 'protein')

write_tsv(orthos, path = 'antismash_cluster_analysis/ETP/cds_regions_focal_genome.orthogroups.long')

orthos %>% 
  group_by(genome, Orthogroup) %>%
  summarize(length(protein)) %>%
  View()
