source('~/scripts/theme_kirsten.R')
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

munge_gca <- function(x){
  paste(strsplit(x, split = "_")[[1]][c(1,2)], collapse = "_")
}

metadata <- read_csv('tables/metadata.csv') %>%
  mutate(acc = sub("\\..*$", "", acc)) %>%
  select(acc, genus_species, Agriculture)

card_files <- list.files('annotation/card', full.names = TRUE)

card <- lapply(card_files, function(x){
  read_tsv(x) %>%
    rename('Resistance_Mechanism' = `Resistance Mechanism`) %>%
    select(ORF_ID, Resistance_Mechanism) %>%
    mutate(genome = ifelse(grepl('GCA', x), yes = munge_gca(basename(x)),
                           no = sub("\\..*$", "", basename(x)))) %>%
    mutate(genome = sub("\\..*$", "", genome))
}) %>% bind_rows()

write_csv(card, path = 'tables/card_annotation.csv')




card_mat <- card %>% 
  group_by(genome, Resistance_Mechanism) %>%
  select(ORF_ID, Resistance_Mechanism) %>%
  summarize(n_genes = length(unique(ORF_ID))) %>%
  filter(n_genes >= 10) %>%
  ungroup() %>%
  spread(Resistance_Mechanism, n_genes) %>% as.data.frame()



rownames(card_mat) <- card_mat$genome
card_mat$genome <- NULL
card_mat <- as.matrix(card_mat)
#colors <- brewer.pal(11, 'Blues')
colors <- colorRampPalette(brewer.pal(n = 7, name ="Blues"))(11)
card_summary %>% filter(n_genes >= 10) %>% .$n_genes %>% summary()
breaks <- seq(min(card_summary$n_genes), max(card_summary$n_genes), by = 30)


pheatmap(card_mat[metadata$acc,],
         breaks = breaks,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 20,
         cellheight = 20,
         color = colors, filename = 'plots/card_heat.pdf')



  
  
  
  
  
  
  
  
  
  
  
  
  

