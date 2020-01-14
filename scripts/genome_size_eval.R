## Command to make the genome scope summary files parsable
## for i in `ls */summary.txt`; do cat $i | sed -e "s/,//g" | sed -e "s/ //" | sed -e "s/ //" | sed -e "s/bp//g" | sed -e "s/%//g">${i%.txt}_fix.txt; done

library('tidyverse')
library('rjson')



jsons <- list.files('kmer_analysis/', pattern = '.dist_analysis.json', full.names = TRUE) 

x <- jsons[1]
kat_gs <- lapply(jsons, function(x){
  json_in <- rjson::fromJSON(file = x)
  json_in$coverage$est_genome_size
})

names(kat_gs) <- jsons %>% basename() %>% sub(".histo.dist_analysis.json", "", .)

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





parent_dir <- 'kmer_analysis'
gs_files <- list.files(parent_dir, 
                       recursive = TRUE,
                       pattern = 'summary_fix.txt')

#x <- gs_files[40]
gs_summary <- lapply(gs_files, function(x){
  read_table2(file.path(parent_dir, x), 
             skip = 6) %>%
    mutate(genome = sub(".genomescope", "", dirname(x))) %>%
    select(-X4) ## strip the empty column<
}) %>% bind_rows() %>%
  mutate(mean = (min+max)/2)
 

# gs_summary %>%
#   select(-min, -max) %>%
#   filter(property %in% c('GenomeUniqueLength', 'GenomeRepeatLength')) %>%
#   ggplot(aes(y = mean, x = genome, fill = property)) +
#     geom_col() + 
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = -45, hjust = 0))





gs_summary %>% filter(property == 'ModelFit') %>% View()

data <- gs_summary %>%
  left_join(., read_coverage, by = 'genome') %>%
  mutate(coverage = bases_total/max) %>%
  group_by(genome) %>%
  mutate(kat_est_gs = kat_gs[[unique(genome)]]) %>%
  ungroup()


data %>%
  select(-bases_total) %>%
  write_csv(., path = 'tables/genome_stats_kmer.csv')

data %>% 
  filter(property == 'GenomeHaploidLength') %>%
  select(genome, min, max, kat_est_gs, coverage) %>%
  write_csv(., path = 'tables/genome_scope_check.csv')


