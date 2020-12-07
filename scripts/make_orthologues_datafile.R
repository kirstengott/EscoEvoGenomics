orthog      <- read_tsv('tables/Orthogroups.tsv') %>% 
  select(-GCA_000225605.1_CmilitarisCM01_v01_protein)

orthog_long <- orthog %>% gather(genome, genes, -Orthogroup) %>%
  mutate(genes = strsplit(genes, split = " ")) %>%
  unnest_longer(genes) %>%
  mutate(genes= sub(",", "", genes)) %>%
  mutate(acc = sub(".all.maker.proteins", '', genome)) %>%
  mutate(acc = sub("\\..*$", "", genome))
  
dir.create('rdata')
save(orthog_long, file = 'rdata/orthologues_long.rdata', compress = 'xz')


annot <- read_tsv('annotation/all_annotations.txt', col_names = c('acc', 'gene', 'tool', 'annot'))


orthog_long_annot <- orthog_long %>%
  rename('gene' = genes) %>%
  left_join(annot, ., by = c('acc', 'gene'))

save(orthog_long_annot, file = 'rdata/orthologues_long_annot.rdata', compress = 'xz')



functional_network <- orthog_long %>%
  rename('gene' = genes) %>%
  left_join(annot, ., by = c('acc', 'gene')) %>% 
  filter(tool %in% c('Pfam')) %>%
  select(-tool, -genome) %>%
  mutate(pfam = strsplit(annot, split = ",")) %>%
  unnest(pfam) %>%
  select(-annot, -gene) %>%
  distinct() %>%
  group_by(acc, pfam) %>%
  mutate(num_ortho = length(unique(Orthogroup))) %>%
  group_by(pfam) %>%
  mutate(num_genomes = length(unique(acc))) %>%
  ungroup()

functional_network <- functional_network %>%
  rename('annot' = pfam) %>% distinct()

save(functional_network, file = 'rdata/functional_network.rdata', compress = 'xz')

gene_network <- orthog_long %>%
  group_by(Orthogroup, genome) %>%
  summarize(num_genes = length(genes))

save(gene_network, file = 'rdata/gene_network.rdata', compress = 'xz')




## get the lengths of the protein annotations
fais <- list.files('annotation/proteins/', pattern = 'fai', full.names = TRUE)
aa_lengths <- lapply(fais, function(x){
  read_tsv(x, col_names = c('gene', 'length', NA, NA, NA)) %>%
    dplyr::select(-contains('X'))
}) %>% bind_rows()



groups = list(c('Higher', 'Leafcutter', 'Lower', 'AllOutgroup'),
              c('Coral', 'Higher', 'Leafcutter', 'AllOutgroup'),
              c('Lower','AllOutgroup'),
              c('Coral', 'Higher', 'Lower', 'AllOutgroup'),
              c('Higher', 'Leafcutter', 'AllOutgroup'),
              c('Coral', 'Lower', 'AllOutgroup'),
              c('Higher', 'Leafcutter'),
              c('Leafcutter'),
              c('Coral', 'Leafcutter', 'Lower', 'AllOutgroup'),
              c('Coral', 'Higher', 'Leafcutter', 'Lower'),
              c('Higher', 'Lower', 'AllOutgroup'),
              c('Leafcutter', 'Lower', 'AllOutgroup'),
              c('Lower'),
              c('Leafcutter', 'AllOutgroup'),
              c('Leafcutter', 'Lower', 'Higher'),
              c('Higher', 'AllOutgroup'),
              c('Coral', 'AllOutgroup'),
              c('Higher'),
              c('Coral', 'Leafcutter', 'Higher'),
              c('Coral', 'Leafcutter', 'AllOutgroup'),
              c('Coral', 'Higher', 'AllOutgroup'),
              c('Lower', 'Higher'),
              c('Lower', 'Leafcutter'),
              c('Coral'),
              c('Coral', 'Lower', 'Higher'),
              c("Coral", 'Lower'),
              c("Coral", 'Leafcutter'),
              c("Coral", 'Lower', 'Leafcutter'),
              c("Coral", 'Higher', 'Leafcutter')
)

groups_sub = list(c('Higher', 'Leafcutter', 'Lower', 'AllOutgroup'),
                  c('Lower','AllOutgroup'),
                  c('Coral', 'Higher', 'Lower', 'AllOutgroup'),
                  c('Coral', 'Higher', 'Leafcutter', 'AllOutgroup'),
                  c('Coral', 'Leafcutter', 'Lower', 'AllOutgroup'),
                  c('Coral', 'Lower', 'AllOutgroup'),
                  c('Higher', 'Leafcutter', 'AllOutgroup'),
                  c('Leafcutter'),
                  c('Lower'),
                  c('Higher', 'Leafcutter'),
                  c('Higher', 'Lower', 'AllOutgroup'),
                  c('Coral', 'Higher', 'Leafcutter', 'Lower'),
                  c('Leafcutter', 'Lower', 'AllOutgroup'))

groups <- lapply(groups, function(x){
  paste(sort(x), collapse = ",")
}) %>% unlist()



input = listInput[c(ag, 'AllOutgroup')]
elements <- unique(unlist(input))


data <- unlist(lapply(input, function(x) {
  x <- as.vector(match(elements, x))
}))

data[is.na(data)] <- as.integer(0)
data[data != 0] <- as.integer(1)
data <- data.frame(matrix(data, ncol = length(input), byrow = F))
data <- data[which(rowSums(data) != 0), ]
names(data) <- paste0("pa_", names(input))
rownames(data) <- elements
data$Orthogroup <- elements




annot_add <- annot %>%
  left_join(., ortho_full, by = 'gene') %>%
  left_join(., aa_lengths, by = 'gene') %>%
  left_join(data, ., by = c('Orthogroup')) %>%
  filter(annot != '-',  
         annot != 'N') %>%
  distinct()
  
  

#save.image(file = 'rdata/OrthologousGeneComparisonData.RData')
