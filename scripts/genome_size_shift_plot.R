library(tidyverse)

raw     <- "jellyfish/kmer_out/kmer_raw_out/LGSR00000000.histo" %>% 
  read_delim(. , delim = " ", comment = "#", col_names = c('occurance', 'count')) %>%
  mutate(type = 'raw')

trimmed <- "jellyfish/kmer_out/kmer_trimmed_out/LGSR00000000.histo" %>% 
  read_delim(. , delim = " ", comment = "#", col_names = c('occurance', 'count')) %>%
  mutate(type = 'trimmed') 

mito    <- "jellyfish/kmer_out/kmer_mito_filtered_out/LGSR00000000.histo" %>% 
  read_delim(. , delim = " ", comment = "#", col_names = c('occurance', 'count')) %>%
  mutate(type = 'mito_filtered')

all_data <- bind_rows(raw, trimmed, mito)
