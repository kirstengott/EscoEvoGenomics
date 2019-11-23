## Command to make the genome scope summary files parsable
## for i in `ls */summary.txt`; do cat $i | sed -e "s/,//g" | sed -e "s/ //" | sed -e "s/ //" | sed -e "s/bp//g" | sed -e "s/%//g">${i%.txt}_fix.txt; done

library('tidyverse')


read_stats <- read_delim('tables/fastp_stats.txt', 
                         col_names = c('file', 'stat', 'value'), 
                         delim = ' ') %>%
  group_by(file) %>%
  mutate(genome = strsplit(file, "_")[[1]][1]) %>%
  ungroup()

read_coverage <- read_stats %>% 
  filter(stat == 'total_bases') %>% 
  group_by(genome) %>%
  summarize(bases_total = sum(value))





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
  mutate(mean = (min+max)/2)
 

gs_summary %>%
  filter(property %in% c('GenomeHaploidLength', 'GenomeRepeatLength')) %>%
  ggplot(aes(y = mean, x = genome, fill = property)) +
    geom_col()





gs_summary %>% filter(property == 'ModelFit') %>% View()

data <- gs_summary %>%
  filter(property == 'GenomeHaploidLength') %>%
  left_join(., read_coverage, by = 'genome') %>%
  mutate(coverage = bases_total/max)

data %>% filter(coverage >= 25) %>% 
  select(genome, min, max, coverage) %>%
  write_csv(., path = 'tables/genome_scope_check.csv')
