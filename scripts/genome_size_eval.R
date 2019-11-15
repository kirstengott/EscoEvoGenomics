library('tidyverse')


read_stats <- read_delim('tables/fastp_stats.txt', col_names = c('genome', 'stat', 'value'), delim = ' ')

read_stats <- read_stats %>% filter(stat == 'total_bases') %>% 
  group_by(genome) %>%
  summarize(bases_total = sum(value))


gs_files <- list.files('tables/jellyfish_histograms/genomescope/', recursive = TRUE,
                       pattern = 'summary_fix.txt')

gs_summary <- lapply(gs_files, function(x){
  read_table2(paste0('tables/jellyfish_histograms/genomescope/', x), 
             skip = 3) %>%
    mutate(genome = dirname(x)) %>%
    select(-X4)
}) %>% bind_rows() %>%
  filter(property == 'GenomeHaploidLength')
 
gs_summary %>% filter(property == 'ModelFit') %>% View()

data <- left_join(gs_summary, read_stats, by = 'genome') %>%
  mutate(min = as.numeric(min)) %>%
  mutate(max = as.numeric(max)) %>%
  mutate(coverage = bases_total/max)

data %>% filter(coverage >= 25) %>% select(genome, min, max, coverage)
