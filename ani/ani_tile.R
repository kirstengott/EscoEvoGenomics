library(tidyverse)
library(pheatmap)
source('~/scripts/theme_kirsten.R')
ani <- read_tsv('ani.txt', col_names = c('r', 'q', 'ani', 'rl', 'ql'), col_types = 'ccnnn') %>%
  mutate(r = basename(r)) %>%
  mutate(q = basename(q))


# ani <- read_tsv('~/Desktop/fastani_out.txt', col_names = c('r', 'q', 'ani', 'rl', 'ql'), col_types = 'ccnnn') %>%
#      mutate(r = sub('.all.maker.transcripts.fasta', '', basename(r))) %>%
#      mutate(q = sub('.all.maker.transcripts.fasta', '', basename(q)))

shorten_genus_species <- function(x){
  spl = strsplit(x, split = " ")
  paste0(strsplit(spl[[1]][1], split = "")[[1]][1], '. ', spl[[1]][2])
}

metadata <- read_csv('metadata.csv') %>%
  mutate(acc_red = sub("\\..*$", "", acc)) %>%
  group_by(genus_species) %>%
  mutate(label = ifelse(is.na(Agriculture), 
                            yes = shorten_genus_species(genus_species),
                            no = genus_species)) %>%
  ungroup()

unique_genomes <- rev(metadata$genome_id)

me_me <- lapply(unique_genomes, function(x){
  c(x, x)
})


#genome_metadata <- read_tsv('qc/genome_metadata.txt')


combinations <- c(combn(x = unique_genomes, 2, simplify = FALSE), me_me)



ani_df <- lapply(combinations, function(t){
  t_matches <- c(paste0(t[1], "_", t[2]), paste0(t[2], "_", t[1]))
  ani %>% mutate(match_id = paste0(r, "_", q)) %>%
    filter(match_id %in% t_matches) %>%
    mutate(ave_ani = mean(ani)) %>%
    mutate(genome1 = t[1]) %>%
    mutate(genome2 = t[2]) %>%
    select(ave_ani, genome1, genome2) %>%
    distinct()
 }) %>%
  bind_rows() %>%
  group_by(genome1) %>%
  mutate(g1 = metadata[which(metadata$genome_id == genome1), 'label', drop = TRUE]) %>%
  group_by(genome2) %>%
  mutate(g2 = metadata[which(metadata$genome_id == genome2), 'label', drop = TRUE]) %>%
  ungroup() %>%
  filter(g1 != g2)  %>%
  mutate(ave_ani_sig = signif(x = ave_ani, digits = 2))





ani_df$g1 <- factor(ani_df$g1, levels = rev(metadata$label))
ani_df$g2 <- factor(ani_df$g2, levels = rev(metadata$label))


ani_df_high <- ani_df %>% filter(ave_ani >= 95)

ggplot(ani_df, aes(x = g1, y = g2, fill = ave_ani)) +
  geom_tile() +
  geom_tile(data = ani_df_high, colour = 'yellow', size = 0.9) +
  #geom_text(data = ani_df_high, aes(label = ave_ani_sig)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave('ani_tile.pdf', dpi = 320, width = 8, height = 10, units = 'in')


ani_df %>% filter(g2 == 'ICBG2049', g1 == 'ICBG736') 

## clades and members with high identity

## 1
## ICBG742, NIGB, LGSR, NQYS

##2
## NQYR, NIGC, NQYQ

## 3
## 1096, 733, 751


## ragoo the genomes of anything that is >99% identical
filter(ani_df, ave_ani_sig >= 98)


#aln_lens <- read_tsv('qc/ANIvis.tsv')
#seq_lens <- read_tsv('qc/fasta_lengths.txt', col_names = c('Reference', 'Reference_Length'))


#all_len <- left_join(aln_lens, seq_lens, by = 'Reference')


## look for things that align over a large length of sequence
#all_len %>% mutate(perc_aligned = (QueryAligned/Reference_Length)*100) %>%
#  filter(perc_aligned > 70)









