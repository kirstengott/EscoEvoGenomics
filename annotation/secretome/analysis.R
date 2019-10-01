library(tidyverse)

shorten_genus_species <- function(x){
  spl = strsplit(x, split = " ")
  paste0(strsplit(spl[[1]][1], split = "")[[1]][1], '. ', spl[[1]][2])
}

## lipases
led_f <- list.files('LED', full.names = TRUE)

x <- led_f[1]



paste(strsplit(basename(x), split = "_")[[1]][c(1,2)], collapse = "_")

led <- lapply(led_f, function(x){
  read_tsv(x, col_names = c('gene', 'lipase', 'eval', 'p_id', 'len')) %>%
    mutate(genome_f = basename(x)) %>%
    mutate(genome = paste(strsplit(genome_f, split = "_")[[1]][c(1,2)], collapse = "_")) %>%
    mutate(genome = sub("\\..*$", "", genome))
}) %>% bind_rows() %>%
  mutate(lipase = 1) %>%
  select(gene, lipase, genome)


led %>% group_by(genome_f) %>%
  summarize(length(unique(lipase)),
            length(unique(gene))) %>%
  View()

## peptidase

merops_f <- list.files('merops', full.names = TRUE)

merops <- lapply(merops_f, function(x){
  read_tsv(x, col_names = c('gene', 'peptidase', 'eval', 'p_id', 'len')) %>%
    mutate(genome_f = basename(x))    %>%
    mutate(genome = paste(strsplit(genome_f, split = "_")[[1]][c(1,2)], collapse = "_")) %>%
    mutate(genome = sub("\\..*$", "", genome))
}) %>% bind_rows() %>%
  mutate(protease = 1) %>%
  select(gene, protease, genome)





## cazymes
fam_db <- read_tsv('CAZyDB.07312018.fam-activities.txt', comment = '#', col_names = c('cazy_base', 'cazy_description'))
cazy_files <- list.files('cazy', full.names = TRUE)

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



interpro_f <- list.files('interproscan/', full.names = TRUE)

interpro <- lapply(interpro_f, function(x){
  read_tsv(x) %>%
    mutate(genome_f = basename(x))    %>%
    mutate(genome = paste(strsplit(genome_f, split = "_")[[1]][c(1,2)], collapse = "_")) %>%
    mutate(genome = sub("\\..*$", "", genome))
}) %>% bind_rows()

## filter IDs that are likely secreted proteins
interpro_filtered <- interpro %>%
  select(genome, Gene, SignalP_EUK, TMHMM) %>%
  filter(!is.na(TMHMM) | !is.na(SignalP_EUK)) %>%
  gather(tool, value, -genome, -Gene) %>%
  filter(!is.na(value)) %>%
  group_by(genome, Gene, value) %>%
  mutate(count = length(value)) %>%
  select(-tool) %>%
  spread(value, count) %>%
  rename('gene' = Gene) %>% .$gene

targetp_f <- list.files('targetp', full.names = TRUE)

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
  distinct() %>%
  select(gene, genome, cazy, cazy_base) %>%
  spread(key = cazy_base, value = cazy)

#cazy_spread[is.na(cazy_spread)] <- 0


all_data <- targetp %>%
  filter(gene %in% interpro_filtered) %>%
  left_join(., cazy_spread, by = c('genome', 'gene')) %>%
  left_join(., merops, by = c('genome', 'gene')) %>%
  left_join(., led, by = c('genome', 'gene')) %>%
  gather(type, value, -gene, -genome) %>%
  mutate(value = ifelse(is.na(value),
                        yes = 0,
                        no = value)) %>%
  group_by(genome, gene) %>%
  mutate(any_activity = sum(value)) %>%
  ungroup() %>%
  filter(any_activity > 0) %>%
  select(-any_activity) %>%
  group_by(genome, type) %>%
  summarize(total = sum(value)) %>%
  ungroup()


## make sure the min is 1, only keep genes with a secretome activity
#rowSums(select(all_data, -gene, -genome)) %>% summary()

# all_data %>% group_by(genome) %>%
#   summarize(length(unique(gene))) %>%
#   View()


metadata <- read_csv('../../Figures/bigscape/metadata.csv') %>%
  mutate(acc_red = sub("\\..*$", "", acc)) %>%
  group_by(genus_species) %>%
  mutate(pca_label = ifelse(is.na(Agriculture),
                            yes = shorten_genus_species(genus_species),
                            no = genus_species)) %>%
  ungroup() %>%
  mutate(acc = sub("\\..*$", "", acc))







## start code for the pca
pca <- all_data %>% spread(key = genome, value = total) %>% data.frame(.)
rownames(pca) <- pca$type
pca$type <- NULL



pca[is.na(pca)] <- 0

pca <- pca[rowSums(pca) > 0, ]

pca_data <- prcomp(pca)


pca_rot <- data.frame(pca_data$rotation,
                      acc = rownames(pca_data$rotation)
) %>%
  left_join(., metadata, by = 'acc') %>%
  mutate(Agriculture = ifelse(is.na(Agriculture),
                              yes = 'Generalist',
                              no = Agriculture)) %>%
  rename('Agriculture' = Agriculture)

colours <- c("Generalist" = "#808080",
             "Coral" = "#CE3DD0",
             "Higher" = "#2D71F6",
             "Lower" = "#FFFEAB",
             "Leafcutter" = "#377D22")

geom_text_data <- pca_rot %>%
  filter(genus != 'Escovopsis')

ggplot(data = pca_rot, aes(x = PC1, y = PC2, label = pca_label)) +
  geom_point(shape = 21, size = 3, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  #geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 5)
  theme_classic() +
  geom_text(data = geom_text_data, hjust = 1, vjust = 0, size = 3) +
  scale_fill_manual(values = colours) +
  scale_x_continuous(expand = c(0, 0.2)) #+
  ggsave('secretome_pca_12.pdf')

ggplot(data = pca_rot, aes(x = PC2, y = PC3, label = pca_label)) +
    geom_point(shape = 21, size = 3, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
    #geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 5)
    theme_classic() +
    geom_text(data = geom_text_data, hjust = 1, vjust = 0, size = 3) +
    scale_fill_manual(values = colours) +
    scale_x_continuous(expand = c(0, 0.2)) +
  ggsave('secretome_pca_23.pdf')
## PCA without trichodermat
keep_accessions <- metadata %>% filter(genus == "Escovopsis") %>% .$acc %>% unique()

## start code for the pca
pca <- all_data %>%
  filter(genome %in% keep_accessions) %>%
  spread(key = genome, value = total) %>% data.frame(.)
rownames(pca) <- pca$type
pca$type <- NULL

pca[is.na(pca)] <- 0

pca<- pca[rowSums(pca) > 0, ]

pca_data <- prcomp(pca)
pca_rot <- data.frame(pca_data$rotation,
                      acc = rownames(pca_data$rotation)
) %>%
  left_join(., metadata, by = 'acc') %>%
  mutate(Agriculture = ifelse(is.na(Agriculture),
                              yes = 'Generalist',
                              no = Agriculture)) %>%
  rename('Agriculture' = Agriculture)

geom_text_data <- pca_rot %>%
  filter(genus != 'Escovopsis')



ggplot(data = pca_rot, aes(x = PC1, y = PC2, label = pca_label)) +
  geom_point(shape = 21, size = 4, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  #geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 5)
  theme_classic() +
  geom_text(data = geom_text_data, hjust = 1, vjust = 0, size = 3) +
  scale_fill_manual(values = colours) +
  scale_x_continuous(expand = c(0, 0.2)) +
  ggsave('secretome_noout_pca_12.pdf')


ggplot(data = pca_rot, aes(x = PC2, y = PC3, label = pca_label)) +
  geom_point(shape = 21, size = 4, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  #geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 5)
  theme_classic() +
  geom_text(data = geom_text_data, hjust = 1, vjust = 0, size = 3) +
  scale_fill_manual(values = colours) +
  scale_x_continuous(expand = c(0, 0.2)) +
  ggsave('secretome_noout_pca_23.pdf')

## all data summary

cazy_summary <- cazy %>% group_by(genome, gene) %>%
  summarize(cazyme = sum(cazy)) %>%
  distinct()

all_data_summary <- targetp %>%
  filter(gene %in% interpro_filtered) %>%
  left_join(., cazy_summary, by = c('genome', 'gene')) %>%
  left_join(., merops, by = c('genome', 'gene')) %>%
  left_join(., led, by = c('genome', 'gene')) %>%
  gather(type, value, -gene, -genome) %>%
  mutate(value = ifelse(is.na(value),
                        yes = 0,
                        no = value)) %>%
  group_by(genome, gene) %>%
  mutate(any_activity = sum(value)) %>%
  ungroup() %>%
  filter(any_activity > 0) %>%
  select(-any_activity) %>%
  group_by(genome, type) %>%
  summarize(total = sum(value)) %>%
  ungroup() %>%
  rename('acc' = genome) %>%
  left_join(., metadata, by = 'acc')

all_data_summary$pca_label <- factor(all_data_summary$pca_label,
                                     levels = rev(metadata$pca_label))


colors <- c('cazyme' = '#46ACC8', 'lipase' = '#E2D200', 'protease' = '#B40F20')

## overall numbers
ggplot(all_data_summary, aes(x = pca_label, y = total, fill = type)) +
  geom_col() +
  labs(x = '', y = '') +
  theme_classic() +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave('secretome_summary.pdf')

## just look at cazy num genes v num domains plot
ag_colours <- c("Generalist" = "#808080",
             "Coral" = "#CE3DD0",
             "Higher" = "#2D71F6",
             "Lower" = "#FFFEAB",
             "Leafcutter" = "#377D22")

cazy_sum <- cazy %>% group_by(genome) %>%
  summarize(num_genes = length(unique(gene)),
            num_domain = length(unique(cazy_base))) %>%
  rename('acc' = genome) %>%
  left_join(., metadata, by = 'acc') %>%
  mutate(Agriculture =ifelse(is.na(Agriculture), yes = 'Generalist', no = Agriculture))

cazy_sum %>%
  ggplot(., aes(x = num_genes, y = num_domain)) +
  geom_point(shape = 21, size = 4, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  scale_fill_manual(values = ag_colours) +
  theme_classic() +
  ggsave('cazy_domain_diversity_all.pdf')


text_sub <- cazy_sum %>% filter(Agriculture != 'Leafcutter', genus == 'Escovopsis')

cazy_sum %>%
  filter(genus == "Escovopsis") %>%
  ggplot(., aes(x = num_genes, y = num_domain, label = pca_label)) +
  geom_point(shape = 21, size = 4, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  geom_text(data = text_sub, hjust = 0, vjust = 0, size = 3) +
  scale_fill_manual(values = ag_colours) +
  scale_x_continuous(expand = c(0.2, 0)) +
  theme_classic() +
  ggsave('cazy_domain_diversity.pdf')


