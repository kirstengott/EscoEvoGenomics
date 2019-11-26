source('~/scripts/theme_kirsten.R')
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(vegan)
library(ggord)

munge_gca <- function(x){
  paste(strsplit(x, split = "_")[[1]][c(1,2)], collapse = "_")
}

metadata <- read_csv('tables/metadata.csv') %>%
  mutate(acc = sub("\\..*$", "", acc)) %>%
  select(acc, genus_species, Agriculture, Host) %>%
  data.frame()
rownames(metadata) <- metadata$acc


card_files <- list.files('annotation/card', full.names = TRUE)

card <- lapply(card_files, function(x){
  read_tsv(x) %>%
    rename('Resistance_Mechanism' = `Resistance Mechanism`) %>%
    select(ORF_ID, Resistance_Mechanism, ID) %>%
    mutate(genome = ifelse(grepl('GCA', x), yes = munge_gca(basename(x)),
                           no = sub("\\..*$", "", basename(x)))) %>%
    mutate(genome = sub("\\..*$", "", genome))
}) %>% bind_rows()

write_csv(card, path = 'tables/card_annotation.csv')

## summarized heatmap per resistance mechanism
card_mat <- card %>% 
  select(-ID) %>%
  group_by(genome, Resistance_Mechanism) %>%
  summarize(n_genes = length(unique(ORF_ID))) %>%
  filter(n_genes >= 10) %>%
  ungroup() %>%
  spread(Resistance_Mechanism, n_genes) %>% as.data.frame()



rownames(card_mat) <- card_mat$genome
card_mat$genome <- NULL
card_mat <- as.matrix(card_mat)
#colors <- brewer.pal(11, 'Blues')
breaks <- quantile(card_mat)
colors <- colorRampPalette(brewer.pal(n = 7, name ="Blues"))(length(breaks))


pheatmap(card_mat[metadata$acc,],
         breaks = breaks,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 20,
         cellheight = 20,
         color = colors, filename = 'plots/card_heat.pdf')


## full ortholog heatmap

card_mat <- card %>%
  group_by(genome, ID) %>%
  summarize(n_genes = length(unique(ORF_ID))) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(tot_genes = sum(n_genes)) %>%
  mutate(tot_genomes = length(unique(genome))) %>%
  ungroup() %>%
  filter(tot_genes >= 10, tot_genomes > 1) %>% ## must have at least ten to look at)
  select(-tot_genes, -tot_genomes) %>%
  spread(ID, n_genes) %>% 
  column_to_rownames(var = 'genome') %>%
  as.matrix()

annotation_row = data.frame(
  "Agriculture" = metadata[rownames(card_mat), 'Agriculture'],
  "Host" = metadata[rownames(card_mat), 'Host'],
  row.names = rownames(card_mat),
  stringsAsFactors = FALSE
)  %>%
  replace_na(list(Agriculture = 'Outgroup', Host = 'Outgroup'))


annotation_col <- filter(card, ID %in% colnames(card_mat)) %>%
  select(Resistance_Mechanism, ID) %>%
  distinct() %>%
  group_by(ID) %>%
  mutate(Resistance = paste(Resistance_Mechanism, collapse = ";")) %>%
  ungroup() %>%
  select(-Resistance_Mechanism) %>%
  distinct() %>%
  column_to_rownames(var = 'ID')




card_mat[is.na(card_mat)] <- 0
breaks <- unique(quantile(card_mat))
colors <- colorRampPalette(brewer.pal(n = 7, name ="Blues"))(length(breaks))


ann_colors = list(
  #'BGC_type' = c(NRPS = '#1B9E77', PKSI = '#D95F02', Terpene = '#7570B3', Others = "#E7298A", PKS.NRP = "#66A61E"),
  'Agriculture' = c(Leafcutter = 'green4', 
                    Lower = 'yellow', 
                    Coral = 'magenta3', 
                    Higher = 'blue3', 
                    Outgroup = 'white')
)



pheatmap(card_mat[metadata$acc,],
         breaks = breaks,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         cellwidth = 20,
         cellheight = 20,
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = FALSE,
         color = colors, filename = 'plots/card_heat_all_orf.pdf')





pca_data <- prcomp(card_mat)

eigs <- pca_data$sdev^2
proportion = (eigs/sum(eigs))*100
head(proportion)
cumulative = cumsum(eigs)/sum(eigs)
screeplot(pca_data)


n_k <- length(which(proportion >= 10))

if(n_k < 2){
  n_k = 2
}

## nmds
## bray curtis distance is more resilient to nulls
example_NMDS=metaMDS(card_mat, k = n_k) # The number of reduced dimensions ## components with >10% variance explained
stressplot(example_NMDS)

treat <- metadata[which(metadata$acc %in%rownames(card_mat)), 'Agriculture', drop = TRUE]
treat <- replace_na(treat, replace = 'Outgroup') ## corresponds to rows/communities (genomes)

## test for significance
dist <- vegdist(card_mat)
ano_test <- anosim(dist, grouping = treat)
summary(ano_test)
plot(ano_test)

## If found no statistical difference amongst agricultural groups
## 

ggord(example_NMDS,
      grp_in = treat,
      arrow = NULL, ## draw the arrows
      obslab = FALSE,
      txt = FALSE,## labeling the ordination
      poly=FALSE, size=2) + theme_classic() +
  labs(subtitle = paste('Pval:', ano_test$signif, ", R:", round(ano_test$statistic, digits = 2))) +
  ggsave('plots/card_ord_all.pdf')

  
  
  
  
  
  
  
  
  
  

