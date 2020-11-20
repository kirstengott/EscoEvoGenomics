orthog      <- read_tsv('tables/Orthogroups.tsv')
orthog_long <- orthog %>% gather(genome, genes, -Orthogroup) %>%
  mutate(genes = strsplit(genes, split = " ")) %>%
  unnest_longer(genes) %>%
  mutate(genes = sub(",", "", genes))
save(orthog_long, file = 'rdata/orthologues_long.rdata', compress = 'xz')
