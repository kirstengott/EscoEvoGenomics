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
  select(acc, genus_species, Agriculture, contains('cazyme_groups')) %>%
  mutate(genus_species = sub(" ", "_", genus_species))


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


## find interesting cazymes by groups
caz_interest <- cazy %>%
  select(-start, -stop, -gene) %>%
  select(-cazy_base, -cazy_description, -Agriculture, -acc, -genus_species, -cazyme_groups1) %>%
  distinct() %>%
  group_by(cazyme_groups2, cazy) %>%
  mutate(n_genes = signif(mean(n_genes), digits =2)) %>%
  ungroup() %>%
  distinct() %>%
  spread(cazyme_groups2, n_genes) %>%
  replace_na(list(Coral1 = 0, FFA = 0, Lower = 0, Outgroup1 = 0, Outgroup2 = 0)) %>%
  gather(treatments, values, -cazy, -Outgroup2) %>%
  mutate(diff = abs(values-Outgroup2)) %>%
  filter(diff >= 5) %>%
  select(-diff) %>%
  spread(treatments, values) %>%
  replace_na(list(Coral1 = 0, FFA = 0, Lower = 0, Outgroup1 = 0, Outgroup2 = 0)) %>%
  data.frame() %>% 
  distinct() %>%
  .$cazy



c_sum <- cazy %>%
  select(-start, -stop, -gene) %>%
  select(-cazy_base, -cazy_description, -Agriculture, -acc, -contains('cazyme_groups')) %>%
  distinct() %>%
  spread(genus_species, n_genes) %>%
  filter(cazy %in% caz_interest) %>%
  data.frame()
   
c_sum[is.na(c_sum)] <- 0

rownames(c_sum) <- c_sum$cazy
c_sum$cazy <- NULL




#c_heat_f <- c_sum[which(abs(rowMeans(c_sum)) >= 3), ]
breaks = c(0, 5, 10, 15, 20, 25, max(c_sum))
colors <- colorRampPalette(brewer.pal(n = 7, name ="Blues"))(length(breaks)-1)

annotation_row = metadata %>% select(genus_species, cazyme_groups2) %>%
  column_to_rownames(var = "genus_species")

pheatmap(as.matrix(t(c_sum)[metadata$genus_species, ]),
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         breaks = breaks,
         border_color = "grey70",
         cellwidth = 15,
         cellheight = 10,
         filename = 'plots/cazy_heat.pdf',
         color = colors,
         display_numbers = TRUE,
         number_format = "%.0f", annotation_row = annotation_row)

c_heat <- cazy %>%
  select(-start, -stop, -gene) %>%
  select(-cazy_base, -cazy_description, -Agriculture, -acc, -contains('cazyme_groups')) %>%
  distinct() %>%
  spread(genus_species, n_genes) %>%
  data.frame()
c_heat[is.na(c_heat)] <- 0

rownames(c_heat) <- c_heat$cazy
c_heat$cazy <- NULL

c_heat <- t(c_heat)

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

treat <- metadata[sapply(labels(c_dist), function(x){grep(x, metadata$genus_species)}), 'cazyme_groups1', drop = TRUE]
treat <- replace_na(treat, replace = 'Outgroup') ## corresponds to rows/communities (genomes)

ano_test <- anosim(c_dist, grouping = treat)
summary(ano_test)
plot(ano_test)

data1 <- example_NMDS$species %>%
  data.frame() %>%
  rownames_to_column(var = 'cazy_base')
data2 <- example_NMDS$points %>% 
  data.frame() %>% 
  rownames_to_column(var = 'genus_species') %>%
  left_join(., metadata, by = 'genus_species')

data1 %>% 
  mutate(interest = ifelse(cazy_base %in% caz_interest, yes = TRUE, no = FALSE)) %>%
  left_join(., fam_db, by = 'cazy_base') %>%
  write_csv('tables/cazy_ord_all.csv')

ggplot(data2, aes(x = MDS1, y = MDS2, color = cazyme_groups1)) + 
  geom_point() + 
#  geom_text(aes(label = genus_species, color = NULL, size = 0.5), nudge_x = 1) + 
  theme_bw() +
  labs(subtitle = paste('Pval:', ano_test$signif, ", R:", signif(ano_test$statistic, digits = 3))) +
  ggsave('plots/cazy_ord_all.pdf')

# ggord(example_NMDS,
#       grp_in = treat,
#       arrow = NULL, ## draw the arrows
#       obslab = FALSE,
#       txt = FALSE,## labeling the ordination
#       poly=FALSE, size=1,
#       ellipse = FALSE) + theme_classic() +
#   labs(subtitle = paste('Pval:', ano_test$signif, ", R:", signif(ano_test$statistic, digits = 3))) +
#   ggsave('plots/cazy_ord_all.pdf')


## only compare leafcutter to the outgroup
treat <- metadata[sapply(labels(c_dist), function(x){grep(x, metadata$genus_species)}), 'cazyme_groups2', drop = TRUE]
treat <- replace_na(treat, replace = 'Outgroup') ## corresponds to rows/communities (genomes)

new_treat_ind <- which(treat %in% c('Coral1', 'Lower', 'FFA'))
new_treat     <- treat[new_treat_ind]
c_new         <- c_heat[new_treat_ind, ]

example_NMDS=metaMDS(c_new, k = n_k) # The number of reduced dimensions ## components with >10% variance explained
stressplot(example_NMDS)

# bgc_dist <- vegdist(c_new)
# ano_test <- anosim(bgc_dist, grouping = new_treat)
# summary(ano_test)
# plot(ano_test)


colors =  c('Lower' = '#FFFEAB',
            'Coral' = "#CE3DD0",
            'Higher' = "#2D71F6",
            'Leafcutter' = "#377D22")

data1 <- example_NMDS$species %>%
  data.frame() %>%
  rownames_to_column(var = 'cazy_base')

data2 <- example_NMDS$points %>% 
  data.frame() %>% 
  rownames_to_column(var = 'genus_species') %>%
  left_join(., metadata, by = 'genus_species')

ggplot(data2, aes(x = MDS1, y = MDS2, color = Agriculture)) + 
  geom_text(data = filter(data2, genus_species %in% c('ICBG712', 'ICBG721')), 
                          aes(label = genus_species, color = NULL), nudge_x = -0.04, hjust = 0) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = colors) +
  labs(subtitle = paste('Pval:', ano_test$signif, ", R:", signif(ano_test$statistic, digits = 3))) +
  ggsave('plots/cazy_ord_sub.pdf')

data1 %>% 
  mutate(interest = ifelse(cazy_base %in% caz_interest, yes = TRUE, no = FALSE)) %>%
  left_join(., fam_db, by = 'cazy_base') %>%
  write_csv('tables/cazy_ord_sub.csv')

# ggord(example_NMDS,
#       grp_in = new_treat,
#       arrow = NULL, ## draw the arrows
#       obslab = FALSE,
#       txt = TRUE,## labeling the ordination
#       poly=FALSE, 
#       size=2, 
#       ellipse = FALSE) + theme_classic() +
#   labs(subtitle = paste('Pval:', ano_test$signif, ", R:", signif(ano_test$statistic, digits = 3))) +
#   ggsave('plots/cazy_ord_sub.pdf')


##TODO
## make a table with comparisons of: compare higher to leafcutter and higher to lower/coral/outgroup

## look at difference between lower/coral/outgroup and leafcutter/higher


data.frame(cazy_base = rownames(c_sum)) %>% 
  left_join(., fam_db, by = 'cazy_base') %>%
   write_csv(., path = 'tables/cazymes_largely_different.csv')


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
