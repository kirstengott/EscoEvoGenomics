keep <- read_csv('tables/genome_scope_checked.csv') %>% 
  filter(histogram_check == 'ok') %>%
  .$genome

parent_dir <- 'kmer_analysis/genomescope/'
gs_files <- list.files(parent_dir, 
                       recursive = TRUE,
                       pattern = 'summary_fix.txt')

gs_summary <- lapply(gs_files, function(x){
  read_table2(paste0(parent_dir, x), 
              skip = 3) %>%
    mutate(genome = dirname(x)) %>%
    select(-X4) ## strip the empty column
}) %>% bind_rows() %>%
  mutate(mean = (min+max)/2) %>% 
  filter(genome %in% keep) %>%
  select(property, genome, mean) %>%
  spread(key = property, value = mean) %>%
  mutate(HaploidLength = GenomeHaploidLength - GenomeRepeatLength) %>%
  gather(property, mean, -genome)


gs_summary %>%
  filter(property %in% c('HaploidLength', 'GenomeRepeatLength')) %>%
  ggplot(aes(y = mean, x = genome, fill = property)) +
  geom_col() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

