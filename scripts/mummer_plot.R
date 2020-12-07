
test <- read_tsv('genome_coverage/mummer/ICBG1054_GCA_004303015.1.coords', col_names = FALSE) %>%
  select(X1, X2, X14) %>%
  rename('start' = X1,
         'end' = X2,
         'chr' = X14) %>% 
  data.frame()
 
test_all <- lapply(rownames(test), function(x){
  data.frame(
    chr = test[x, 'chr'],
    pos_all = seq(from = test[x, 'start'], to = test[x, 'end']), 
    stringsAsFactors = FALSE)
}) %>% bind_rows()


index <- read_tsv('genome_coverage/genome_index/GCA_004303015.1_ASM430301v1_genomic.fna.fai',
                  col_names = c('chr', 'length' , NA, NA)) %>% 
  select(-contains('X')) %>%
  data.frame()


seq_all <- lapply(rownames(index), function(x){
  data.frame(
    chr = index[x, 'chr'],
    pos_all = seq(1, index[x, 'length']), 
             stringsAsFactors = FALSE)
}) %>% bind_rows()



all <- left_join(seq_all, test_all, by = c('chr', 'pos_all'))


test <- read_tsv('genome_coverage/mummer/ICBG1054_GCA_004303015.1.coords', col_names = FALSE) %>%
  select(X1, X2, X14) %>%
  rename('start' = X1,
         'end' = X2,
         'chr' = X14)

ggplot(test, aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5)) + 
  geom_rect() + 
  theme_classic() +
  facet_wrap(~chr, nrow = 1)

