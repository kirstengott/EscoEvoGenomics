library(tidyverse)

dn_ds_f <- read_csv('tables/dn_ds_all.txt', col_names = c('ortho', 'g1', 'g2', 'dn_ds'))
orthog  <- read_tsv('tables/Orthogroups.tsv')

busco_f <-list.files('annotation/ascomybocta_odb10_busco', full.names = TRUE)


genomes <- colnames(orthog) %>% 
  purrr::map(.x = ., function(x){sub(pattern = ".all.maker.proteins", replacement = "", x)}) %>%
  purrr::map(.x = ., function(x){sub(pattern = ".protein.", replacement = "", x)}) %>%
  purrr::map(.x = ., function(x){sub(pattern = ".protein", replacement = "", x)}) %>%
  unlist()
genomes <- genomes[-c(1, 5)]

index <- sapply(genomes, function(x){
  grep(x, busco_f)
})

 busco <- lapply(seq(1, length(index)), function(x){
   genome = names(index[x])
   read_tsv(busco_f[index[x]], skip = 3, col_names = c('busco', 'status',	'sequence', 'score', 	'length')) %>%
     mutate(genome = genome)
}) %>% bind_rows() %>%
   filter(!is.na(sequence))

card_f <- list.files('annotation/card', full.names = TRUE)

index <- sapply(genomes, function(x){
  grep(x, card_f)
})


card <- lapply(seq(1, length(index)), function(x){
  genome = names(index[x])
  read_tsv(card_f[index[[x]]]) %>%
    mutate(genome = genome)
}) %>% bind_rows()


