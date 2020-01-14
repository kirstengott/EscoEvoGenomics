## for CAFE
# The data file may consist of data on one family or up to thousands of families for the specified tree.
# The first line of the data file should contain the extant species’ names (as used in the Newick tree description),
# tab-delimited in no particular order.
# Subsequent lines each correspond to a gene family and contain tab-delimited family sizes for these extant species.
# Columns in the data file whose header does not correspond to any of the names in the Newick tree description,
# which may provide additional information about the gene families, are ignored and simply copied in the output file


cazy_caf <- cazy %>% select(genome, cazy_base, cazy_description, n_genes) %>% distinct() %>%
  #  spread(key = genome, value = n_genes) %>%
  rename('family' = cazy_base,
         'family_description' = cazy_description)
glimpse(cazy_caf)

bigscape_caf <- bigscape %>% mutate(Genome = sub("\\..*$", "", Genome)) %>%
  #spread(key = Genome, value = total) %>%
  rename('family' = Type,
         'genome' = Genome,
         'n_genes' = total) %>%
  mutate(family_description = 'null') #%>% select(family, family_description, everything())

glimpse(bigscape_caf)

card_caf <- card %>% mutate(genome = sub("\\..*$", "", genome)) %>%
  #spread(key = genome, value = n_genes) %>%
  rename('family' = Resistance_Mechanism) %>%
  mutate(family_description = 'null') #%>% select(family, family_description, everything())

glimpse(card_caf)

metadata <- read_tsv('id_metadata.txt') %>%
  mutate(genome = sub("\\..*$", "", acc)) %>%
  select(-acc, -family, -genus, -genome_id, -id)

all_cafe <- bind_rows(cazy_caf, bigscape_caf, card_caf) %>%
  ungroup() %>%
  left_join(., metadata, by = 'genome') %>%
  select(-genome) %>%
  spread(key = genus_species, value = n_genes)


all_cafe[is.na(all_cafe)] <- 0


glimpse(all_cafe)

write_tsv(all_cafe, path = 'cafe_table_input.tsv')


