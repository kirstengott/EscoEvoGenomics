library(tidyverse)
setwd("~/Analysis/tree/cds3")
tax <- read_csv('../taxonomy.csv', col_names = FALSE) %>%
  rename('file' = X4) %>%
  separate(X7, into = c(rep(NA, 12), 'order', 'family' ,'genus', 'genus_species'), sep = ";") %>%
  select(-X8:-X24) %>%
  rename('gbrs_paired_asm' = X1) %>%
  mutate(file = sub("\\..*$", "", file))


refseq <- read_tsv('refseq_metadata.txt') %>%
  mutate(file = sub("\\..*$", "", basename(local_filename))) %>%
  left_join(., tax, by = 'gbrs_paired_asm') %>%
  rename('file' = file.x)

tree <- read_tsv('id_map.txt', col_names = c('id', 'file')) %>%
  mutate(file = sub("\\..*$", "", file)) %>%
  left_join(., refseq, by = 'file') %>%
   mutate(genus = case_when(
     file == 'GCA_000326105' ~ 'saccharomyces',
     file == 'GCA_004303015' ~ 'cladobotryum',
     file == 'GCA_000149205' ~ 'aspergillus',
     file == 'GCA_011799845' ~ 'hypomyces',
     file == 'SPDT00000000' ~ 'hypomyces',
     is.na(genus) ~ 'escovopsis',
     TRUE ~ genus
   ))

glimpse(tree)
nrow(tree)






tree %>% write_csv(., path = 'tree_metadata.csv', col_names = TRUE)


