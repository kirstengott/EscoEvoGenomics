source('~/scripts/theme_kirsten.R')
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(ggfortify)
library(vegan)
library(ggord)

cazy_files <- list.files('annotation/cazy', full.names = TRUE)
metadata <- read_csv('tables/metadata.csv') %>%
  mutate(acc = sub("\\..*$", "", acc)) %>%
  select(acc, genus_species, Agriculture, cazyme_groups)


fam_db <- read_tsv('annotation/CAZyDB.07312018.fam-activities.txt', comment = '#', col_names = c('cazy_base', 'cazy_description'))

cazy <- lapply(cazy_files, function(c){
  read_tsv(c) %>%
    mutate(genome = sub(".txt", "", basename(c)))
}) %>% bind_rows() %>%
  separate(HMMER, into = c('hmmer', 'startstop'), sep = "\\(", remove = FALSE) %>%
  separate(startstop, into = c('start', 'stop'), sep = '-') %>%
  mutate(stop = sub("\\)", "", stop)) %>%
  mutate(HMMER = gsub('\\([0-9-]+\\)', "", HMMER)) %>%
  select(-hmmer) %>%
  rename(gene = `Gene ID`) %>%
  gather(tool, cazy_id, -gene, -start, -stop, -`#ofTools`, -genome) %>%
  group_by(gene) %>%
  mutate(cazy_ids = paste(unique(cazy_id), collapse = ",")) %>%
  select(-tool, -cazy_id, -`#ofTools`) %>%
  unnest(cazy_id = strsplit(cazy_ids, ",")) %>%
  unnest(cazy = strsplit(cazy_id, "\\+")) %>%
  select(gene, start, stop, genome, cazy) %>%
  distinct() %>%
  mutate(cazy_base = sub("_.*$", "", cazy)) %>%
  left_join(., fam_db, by = 'cazy_base') %>%
  filter(cazy != 'N') %>%
  group_by(genome, cazy_base) %>%
  mutate(n_genes = as.numeric(length(unique(gene)))) %>%
  ungroup() %>%
  rename('acc' = genome) %>%
  left_join(., metadata, by = 'acc')


write_csv(cazy, path = 'tables/cazy_annotation.csv')



c_sum <- cazy %>%
  select(-start, -stop, -gene) %>%
  select(-cazy_base, -cazy_description, -Agriculture, -acc, -genus_species) %>%
  distinct() %>%
  group_by(cazyme_groups, cazy) %>%
  mutate(n_genes = mean(n_genes)) %>%
  ungroup() %>%
  distinct() %>%
  spread(cazyme_groups, n_genes)

c_heat <- c_sum  %>%
  data.frame()
c_heat[is.na(c_heat)] <- 0

rownames(c_heat) <- c_heat$cazy
c_heat$cazy <- NULL


c_heat_f <- c_heat[which(rowSums(c_heat) > 20), ]

colors <- colorRampPalette(brewer.pal(n = 7, name ="Blues"))(4)
## think about taking rolling differences ?
# this heatmap isn't very helpful with too much data
pheatmap(as.matrix(c_heat),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border_color = "grey70",
         cellwidth = 10,
         cellheight = 10,
         filename = 'cazy_heat.pdf',
         color = colors)


pca_data <- prcomp(c_heat)

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
example_NMDS=metaMDS(c_heat, k = n_k) # The number of reduced dimensions ## components with >10% variance explained
stressplot(example_NMDS)



## test for significance
c_dist <- vegdist(c_heat)

treat <- metadata[sapply(labels(c_dist), function(x){grep(x, metadata$genus_species)}), 'Agriculture', drop = TRUE]
treat <- replace_na(treat, replace = 'Outgroup') ## corresponds to rows/communities (genomes)

ano_test <- anosim(c_dist, grouping = treat)
summary(ano_test)
plot(ano_test)



ggord(example_NMDS,
      grp_in = treat,
      arrow = NULL, ## draw the arrows
      obslab = TRUE,
      txt = FALSE,## labeling the ordination
      poly=FALSE, size=2,
      ellipse = FALSE) + theme_classic() +
  labs(subtitle = paste('Pval:', ano_test$signif, ", R:", round(ano_test$statistic, digits = 2))) +
  ggsave('plots/cazy_ord_all.pdf')


## only compare leafcutter to the outgroup
new_treat_ind <- which(treat %in% c('Higher', 'Coral', 'Lower', 'Leafcutter'))
new_treat     <- treat[new_treat_ind]
c_new         <- c_heat[new_treat_ind, ]

example_NMDS=metaMDS(c_new, k = n_k) # The number of reduced dimensions ## components with >10% variance explained
stressplot(example_NMDS)
bgc_dist <- vegdist(c_new)
ano_test <- anosim(bgc_dist, grouping = new_treat)
summary(ano_test)
plot(ano_test)


ggord(example_NMDS,
      grp_in = new_treat,
      arrow = NULL, ## draw the arrows
      obslab = TRUE,
      txt = FALSE,## labeling the ordination
      poly=FALSE, 
      size=2, 
      ellipse = FALSE) + theme_classic() +
  labs(subtitle = paste('Pval:', ano_test$signif, ", R:", round(ano_test$statistic, digits = 2))) +
  ggsave('plots/cazy_ord_sub.pdf')


##TODO
## make a table with comparisons of: compare higher to leafcutter and higher to lower/coral/outgroup

## look at difference between lower/coral/outgroup and leafcutter/higher


nmds_table <- data.frame(example_NMDS$species, stringsAsFactors = FALSE)
nmds_table$cazy_base <- rownames(nmds_table)


# nmds_table %>% left_join(., fam_db, by = 'cazy_base') %>%
#   write_tsv(., path = 'cazy_nmds_table.tsv')


#pca_data <- prcomp(c_heat)
#
# metadata_t <- metadata %>% mutate(acc = sub("\\..*$", "", acc))
#
# pca_rot <- data.frame(pca_data$rotation,
#                       acc = rownames(pca_data$rotation)
# ) %>%
#   left_join(., metadata_t, by = 'acc') %>%
#   mutate(Agriculture = ifelse(is.na(Agriculture),
#                               yes = 'Generalist',
#                               no = Agriculture)) %>%
#   rename('Agriculture' = Agriculture)
#
# colours <- c("Generalist" = "#808080",
#              "Coral" = "#CE3DD0",
#              "Higher" = "#2D71F6",
#              "Lower" = "#FFFEAB",
#              "Leafcutter" = "#377D22")
#
# geom_text_data <- pca_rot %>%
#   filter(genus != 'Escovopsis')
#
# ggplot(data = pca_rot, aes(x = PC1, y = PC2)) +
#   geom_point(shape = 21, size = 4, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
#   #geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 5)
#   theme_classic() +
#   #geom_text(data = geom_text_data, hjust = 1, vjust = 0, size = 3) +
#   scale_fill_manual(values = colours) +
#   scale_x_continuous(expand = c(0, 0.2)) +
#   ggsave('cazy_pca_12.pdf')
#
# eigs <- pca_data$sdev^2
# proportion = (eigs/sum(eigs))*100
# cumulative = cumsum(eigs)/sum(eigs)
# ## look at this
# head(pca_data$x)
# ## bigscape PCA
#
# autoplot(pca_data, data = c_heat,
#          loadings = TRUE, loadings.colour = 'blue',
#          loadings.label = TRUE, loadings.label.size = 3) + ggsave('test.pdf')
