library(tidyverse)
library(vegan)

## just look at cazy num genes v num domains plot
ag_colours <- c("Generalist" = "#808080",
                "Coral" = "#CE3DD0",
                "Higher" = "#2D71F6",
                "Lower" = "#FFFEAB",
                "Leafcutter" = "#377D22",
                "Outgroup1" = 'black',
                "Outgroup2" = 'gray')


shorten_genus_species <- function(x){
  spl = strsplit(x, split = " ")
  paste0(strsplit(spl[[1]][1], split = "")[[1]][1], '. ', spl[[1]][2])
}

## lipases
led_f <- list.files('annotation/LED', full.names = TRUE)

x <- led_f[1]

paste(strsplit(basename(x), split = "_")[[1]][c(1,2)], collapse = "_")

led <- lapply(led_f, function(x){
  read_tsv(x, col_names = c('gene', 'lipase', 'eval', 'p_id', 'len')) %>%
    mutate(genome_f = basename(x)) %>%
    mutate(genome = paste(strsplit(genome_f, split = "_")[[1]][c(1,2)], collapse = "_")) %>%
    mutate(genome = sub("\\..*$", "", genome))
}) %>% bind_rows() %>%
  select(gene, lipase, genome)


led %>% group_by(genome_f) %>%
  summarize(length(unique(lipase)),
            length(unique(gene))) %>%
  View()

## peptidase

merops_f <- list.files('annotation/merops', full.names = TRUE)

merops <- lapply(merops_f, function(x){
  read_tsv(x, col_names = c('gene', 'peptidase', 'eval', 'p_id', 'len')) %>%
    mutate(genome_f = basename(x))    %>%
    mutate(genome = paste(strsplit(genome_f, split = "_")[[1]][c(1,2)], collapse = "_")) %>%
    mutate(genome = sub("\\..*$", "", genome)) %>%
    mutate(peptidase = sub("-.*$", "", peptidase))
}) %>% bind_rows() %>%
  select(gene, peptidase, genome)





## cazymes
fam_db <- read_tsv('annotation/CAZyDB.07312018.fam-activities.txt', comment = '#', col_names = c('cazy_base', 'cazy_description'))
cazy_files <- list.files('annotation/cazy', full.names = TRUE)

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
  mutate(cazy = 1) #%>% select(gene, cazy, genome)



# interpro_f <- list.files('annotation/interpro/', full.names = TRUE)
# 
# interpro <- lapply(interpro_f, function(x){
#   read_tsv(x) %>%
#     mutate(genome_f = basename(x))    %>%
#     mutate(genome = paste(strsplit(genome_f, split = "_")[[1]][c(1,2)], collapse = "_")) %>%
#     mutate(genome = sub("\\..*$", "", genome))
# }) %>% bind_rows()
# 
# ## filter IDs that are likely secreted proteins
# interpro_filtered <- interpro %>%
#   select(genome, Gene, SignalP_EUK, TMHMM) %>%
#   filter(!is.na(TMHMM) | !is.na(SignalP_EUK)) %>%
#   gather(tool, value, -genome, -Gene) %>%
#   filter(!is.na(value)) %>%
#   group_by(genome, Gene, value) %>%
#   mutate(count = length(value)) %>%
#   select(-tool) %>%
#   spread(value, count) %>%
#   rename('gene' = Gene) %>%


targetp_f <- list.files('annotation/targetp', full.names = TRUE)

# Prediction of localization, based on the scores above; the possible values are:
#   C	Chloroplast, i.e. the sequence contains cTP, a chloroplast transit peptide;
# M	Mitochondrion, i.e. the sequence contains mTP, a mitochondrial targeting peptide;
# S	Secretory pathway, i.e. the sequence contains SP, a signal peptide;

targetp <- lapply(targetp_f, function(x){
  read_tsv(x) %>%
    mutate(genome_f = basename(x))    %>%
    mutate(genome = paste(strsplit(genome_f, split = "_")[[1]][c(1,2)], collapse = "_")) %>%
    mutate(genome = sub("\\..*$", "", genome))
}) %>% bind_rows() %>%
  filter(Loc == 'S') %>% ## pull out secreted proteins
  rename('gene' = Name) %>%
  filter(Len <= 300) %>%
  select(gene, genome)


#wolfp_f <- list.files('wolfpsort', full.names = TRUE)

cazy_spread <- cazy %>%
  select(gene, genome, cazy_base) %>%
  rename('cazyme' = cazy_base)

#cazy_spread[is.na(cazy_spread)] <- 0
all_annotations <- read_tsv('annotation/all_annotations.txt', col_names = c('genome', 'gene', 'tool', 'value'))

## number of annotations
all_annotations %>% select(genome, gene) %>% distinct() %>% .$genome %>% table()

interpro_filtered <- all_annotations %>% 
  filter(tool %in% c('TMHMM', 'SignalP_EUK')) %>% 
  tidyr::spread(key = tool, value = value) %>% 
  filter(!is.na(SignalP_EUK) | TMHMM == 1) %>%
  .$gene


## I'm looking for extr

## cytoskeleton, cytoplasm, nucleus, mitochondria, vesicles of secretory system, 
## endoplasmic reticulum (ER), Golgi, vacuole, plasma membrane, peroxisome, 
## extracellular space including cell wall.

wolfpsort <- all_annotations %>% filter(tool == 'wolfpsort', grepl('extr', value)) %>%
  select(gene, genome) %>%
  mutate(genome = sub("\\..*$", "", genome))


all_secreted <- bind_rows(targetp, wolfpsort) %>% distinct()

all_enzymes <- cazy_spread %>%
  left_join(., merops, by = c('genome', 'gene')) %>%
  left_join(., led, by = c('genome', 'gene')) %>%
  gather(type, value, -gene, -genome) %>%
  filter(!is.na(value)) %>%
  mutate(count = 1)


secreted <- all_secreted %>%
  filter(gene %in% interpro_filtered) %>%
  left_join(., all_enzymes, by = c('genome', 'gene')) %>%
  filter(!is.na(type))



secreted_nmds <- secreted %>%
  select(genome, value, count) %>%
  group_by(genome, value) %>%
  summarize(count_all = sum(count)) %>%
  spread(key = value, value = count_all) %>%
  column_to_rownames(var  = 'genome') %>%
  data.frame()
secreted_nmds[is.na(secreted_nmds)] <- 0


pca_data <- prcomp(secreted_nmds)
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
example_NMDS=metaMDS(secreted_nmds, k = n_k) # The number of reduced dimensions ## components with >10% variance explained
stressplot(example_NMDS)

## test for significance
c_dist <- vegdist(secreted_nmds)

treat <- metadata[sapply(labels(c_dist), function(x){grep(x, metadata$genome_id)}), 'Agriculture', drop = TRUE]

#treat <- replace_na(treat, replace = 'Outgroup') ## corresponds to rows/communities (genomes)

ano_test <- anosim(c_dist, grouping = treat)
summary(ano_test)
plot(ano_test)

data1 <- example_NMDS$species %>%
  data.frame() %>%
  rownames_to_column(var = 'cazy_base')

data2 <- example_NMDS$points %>% 
  data.frame() %>% 
  rownames_to_column(var = 'acc')

data2 <- metadata[sapply(data2$acc, function(x){grep(x, metadata$genome_id)}),] %>%
bind_cols(., data2)

data2_sub <- data2 %>% filter(genus_species %in% c('Cladobotryum protrusum',
                                                   'Hypomyces perniciosus',
                                                   'ICBG712', 
                                                   'ICBG721'))

ggplot(data2, aes(x = MDS1, y = MDS2)) + 
  geom_point(shape = 21, size = 3, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  #geom_text(data = data1_sub, aes(label = cazy_base, color = NULL)) +
  scale_fill_manual(values = ag_colours) +
  scale_color_manual(values = ag_colours) +
  geom_text(data= data2_sub, aes(label = genus_species, 
                                 color = NULL,
                                 hjust=0, vjust=0)) + 
  theme_bw() +
  ggsave('plots/secretome_ord_all.pdf', width = 7, height = 7)





# all_secreted_summary$genus_species <- factor(all_secreted_summary$genus_species,
#                                              levels = rev(metadata$genus_species))
# 
# 
# colors <- c('cazyme' = '#46ACC8', 'lipase' = '#E2D200', 'peptidase' = '#B40F20')
# 
# ## overall numbers
# ggplot(all_secreted_summary, aes(x = genus_species, y = count, fill = type)) +
#   coord_flip() +
#   geom_col() +
#   labs(x = '', y = '') +
#   theme_classic() +
#   scale_fill_manual(values = colors) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   ggsave('plots/secretome_summary.pdf')



not_secreted <- all_enzymes %>% filter(!gene %in% secreted$gene)

not_secreted_nmds <- not_secreted %>%
  select(genome, value, count) %>%
  group_by(genome, value) %>%
  summarize(count_all = sum(count)) %>%
  spread(key = value, value = count_all) %>%
  column_to_rownames(var  = 'genome') %>%
  data.frame()
not_secreted_nmds[is.na(not_secreted_nmds)] <- 0


pca_data <- prcomp(not_secreted_nmds)
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
example_NMDS=metaMDS(not_secreted_nmds, k = n_k) # The number of reduced dimensions ## components with >10% variance explained
stressplot(example_NMDS)

## test for significance
c_dist <- vegdist(not_secreted_nmds)

treat <- metadata[sapply(labels(c_dist), function(x){grep(x, metadata$genome_id)}), 'Agriculture', drop = TRUE]

#treat <- replace_na(treat, replace = 'Outgroup') ## corresponds to rows/communities (genomes)

ano_test <- anosim(c_dist, grouping = treat)
summary(ano_test)
plot(ano_test)


data2 <- example_NMDS$points %>% 
  data.frame() %>% 
  rownames_to_column(var = 'acc')

data2 <- metadata[sapply(data2$acc, function(x){grep(x, metadata$genome_id)}),] %>%
  bind_cols(., data2)

data2_sub <- data2 %>% filter(genus_species %in% c('Cladobotryum protrusum',
                                                   'Hypomyces perniciosus',
                                                   'ICBG712', 
                                                   'ICBG721'))

ggplot(data2, aes(x = MDS1, y = MDS2)) + 
  geom_point(shape = 21, size = 3, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  #geom_text(data = data1_sub, aes(label = cazy_base, color = NULL)) +
  scale_fill_manual(values = ag_colours) +
  geom_text(data= data2_sub, aes(label = genus_species, 
                                 color = NULL,
                                 hjust=0, vjust=0)) + 
  theme_bw() +
  ggsave('plots/non_secretome_ord_all.pdf', width = 7, heigh = 7)


all_secreted_summary <- secreted %>% 
  group_by(genome, type) %>%
  summarize(count = sum(count))%>%
  mutate(group = 'secreted')

catabolite_summary <- not_secreted %>% 
  group_by(genome, type) %>%
  summarize(count = sum(count)) %>%
  mutate(group = 'non_secreted') %>%
  bind_rows(., all_secreted_summary) %>%
  rename('acc' = genome) %>%
  left_join(., metadata, by = 'acc') 

catabolite_summary$genus_species <- factor(catabolite_summary$genus_species,
                                             levels = rev(metadata$genus_species))

## overall numbers
ggplot(catabolite_summary, aes(x = genus_species, y = count, fill = type)) +
  coord_flip() +
  geom_col() +
  labs(x = '', y = '') +
  facet_grid(~group) +
  theme_classic() +
  scale_fill_manual(values = colors) +
  ggsave('plots/catabolite_summary.pdf', width = 14)



metadata <- read_csv('tables/metadata.csv') %>%
  mutate(acc_red = sub("\\..*$", "", acc)) %>%
  group_by(genus_species) %>%
  mutate(pca_label = ifelse(is.na(Agriculture),
                            yes = shorten_genus_species(genus_species),
                            no = genus_species)) %>%
  ungroup() %>%
  mutate(acc = sub("\\..*$", "", acc))




all_sum <- all_enzymes %>% 
  group_by(genome, type) %>%
  summarize(num_genes = length(unique(gene)),
            num_domain = length(unique(value))) %>%
  rename('acc' = genome) %>%
  left_join(., metadata, by = 'acc') %>%
  mutate(Agriculture =ifelse(is.na(Agriculture), 
                             yes = 'Generalist', 
                             no = Agriculture))

all_sum %>%
  ggplot(., aes(x = num_genes, y = num_domain)) +
  geom_point(shape = 21, size = 3, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  scale_fill_manual(values = ag_colours) +
  facet_grid(cols = vars(type), scales = 'free_x') + 
  scale_y_continuous(breaks = seq(0, 150, by = 10), labels =seq(0, 150, by = 10)) +
  theme_bw() +
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  ) + 
  ggsave('plots/catabolic_domain_diversity_all.pdf', width = 15, height = 8)


## enzymes only in ffa
all_enzymes %>%
  rename('acc' = genome) %>%
  left_join(., metadata, by = 'acc') %>%
  mutate(Agriculture =ifelse(is.na(Agriculture), 
                             yes = 'Generalist', 
                             no = Agriculture)) %>%
  select(cazyme_groups1, value) %>%
  distinct() %>%
  group_by(value) %>%
  mutate(num_groups = length(unique(cazyme_groups1))) %>%
  filter(num_groups == 1 && cazyme_groups1 == 'FFA')
  


