# ls *GCA_004303015.1*filtered.bam | tr '\n' '\t' | sed -e "s/^/chr\tpos\t/" | sed -e 's/_GCA_004303015.1.sorted.filtered.bam//g' >GCA_004303015.1.depth
# samtools depth *GCA_004303015.1*filtered.bam >>GCA_004303015.1.depth



library(tidyverse)

depth <- read_tsv('genome_coverage/depth/SPDT00000000.1.depth.gz')

index <- read_tsv('genome_coverage/genome_index/SPDT00000000.1_genomic.fna.fai',
                  col_names = c('chr', 'length' , NA, NA)) %>% 
  select(-contains('X'))

## millions of reads that aligned
sum(depth$count)/1000000

depth = depth[,which(colSums(depth[, -c(1,2)]) > 0)]

depth_f <- depth %>% gather(genome, count, -chr, -pos) %>%
  mutate(adj_count = ifelse(count < 10, 0, count))




lapply(depth$chr, function(x){
  depth_f %>% filter(chr == x,
                     genome %in% c('751', 
                                   'ICBG1065',
                                   'ICBG2046')) %>% 
    ggplot(., aes(y = adj_count, x = pos)) +
    geom_line() +
    facet_wrap(~genome, ncol = 1) +
    theme_classic()
})











