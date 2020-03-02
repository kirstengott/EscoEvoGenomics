library(tidyverse)

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


# led %>% group_by(genome_f) %>%
#   summarize(length(unique(lipase)),
#             length(unique(gene))) %>%
#   View()

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


viru_f <- list.files('annotation/virulence/parsed', full.names = TRUE, pattern = 'fa')

viru <- lapply(viru_f, function(x){
  read_tsv(x, col_names = c('gene', 'virulence', 'eval', 'p_id', 'len')) %>%
    mutate(genome_f = basename(x))    %>%
    mutate(genome = paste(strsplit(genome_f, split = "_")[[1]][c(1,2)], collapse = "_")) %>%
    mutate(genome = sub("\\..*$", "", genome)) %>%
    mutate(peptidase = sub("-.*$", "", virulence))
}) %>% bind_rows() %>%
  select(gene, virulence, genome)


## cazymes
fam_db <- read_tsv('annotation/CAZyDB.07312018.fam-activities.txt', comment = '#', col_names = c('cazy_base', 'cazy_description'))

cazy_files <- list.files('annotation/cazy/eCAMIout', full.names = TRUE)

cazy <- lapply(cazy_files, function(c){
  read_tsv(c, col_names = c('gene', 'family', 'subfamily')) %>%
    mutate(gene = sub(" .*$", "", gene)) %>%
    mutate(genome = basename(c))    %>%
    mutate(genome = paste(strsplit(genome, split = "_")[[1]][c(1,2)], collapse = "_")) %>%
    mutate(genome = sub("\\..*$", "", genome)) %>%
    mutate(cazyme = gsub('\\|[0-9]\\..*$', '', subfamily)) %>% ## remove counts    
    mutate(cazyme = sub("\\|.*$", '', cazyme)) %>%
    mutate(cazyme = gsub(':[0-9]+', '', cazyme)) %>% ## remove counts
    mutate(ec = stringr::str_extract_all(subfamily, '\\|[0-9\\.:]+\\|')) %>%
    mutate(ec = sapply(ec, paste, collapse='')) %>%
    mutate(ec = gsub(':[0-9]+', '', ec)) %>%
    mutate(ec = gsub('\\|\\|', ",", ec)) %>%
    mutate(ec = gsub('\\|', '', ec))## remove counts
}) %>% bind_rows() %>%
  distinct() %>%
  group_by(genome, cazyme) %>%
  mutate(n_genes = as.numeric(length(unique(gene)))) %>%
  ungroup() %>%
  group_by(genome, gene) %>%
  mutate(n_cazy = as.numeric(length(unique(cazyme)))) %>%
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
  select(gene, genome, cazyme)

#cazy_spread[is.na(cazy_spread)] <- 0
all_annotations <- read_tsv('annotation/all_annotations.txt', col_names = c('genome', 'gene', 'tool', 'value'))

interpro_filtered <- all_annotations %>% 
  filter(tool %in% c('TMHMM', 'SignalP_EUK')) %>% 
  tidyr::spread(key = tool, value = value) %>% 
  filter(!is.na(SignalP_EUK) | TMHMM == 1) %>%
  .$gene

wolfpsort <- all_annotations %>% filter(tool == 'wolfpsort', grepl('extr', value)) %>%
  select(gene, genome) %>%
  mutate(genome = sub("\\..*$", "", genome))


all_secreted <- bind_rows(targetp, wolfpsort) %>% distinct()

all_enzymes <- cazy_spread %>%
  left_join(., merops, by = c('genome', 'gene')) %>%
  left_join(., led, by = c('genome', 'gene')) %>%
  left_join(., viru, by = c('genome', 'gene')) %>%
  gather(type, value, -gene, -genome) %>% 
  filter(!is.na(value)) %>%
  mutate(count = 1)


secreted <- all_secreted %>%
  filter(gene %in% interpro_filtered) %>%
  left_join(., all_enzymes, by = c('genome', 'gene')) %>%
  filter(!is.na(type))

## number of annotations
all_annotations %>% select(genome, gene) %>% distinct() %>% .$genome %>% table()


## I'm looking for extr

## cytoskeleton, cytoplasm, nucleus, mitochondria, vesicles of secretory system, 
## endoplasmic reticulum (ER), Golgi, vacuole, plasma membrane, peroxisome, 
## extracellular space including cell wall.


## make sure the min is 1, only keep genes with a secretome activity
#rowSums(select(all_data, -gene, -genome)) %>% summary()

# all_data %>% group_by(genome) %>%
#   summarize(length(unique(gene))) %>%
#   View()



all_secreted_summary <- secreted %>% 
  group_by(genome, type) %>%
  summarize(count = sum(count))
  
  
metadata <- read_csv('tables/metadata.csv') %>%
  mutate(acc_red = sub("\\..*$", "", acc)) %>%
  group_by(genus_species) %>%
  mutate(pca_label = ifelse(is.na(Agriculture),
                            yes = shorten_genus_species(genus_species),
                            no = genus_species)) %>%
  ungroup() %>%
  mutate(acc = sub("\\..*$", "", acc))



colors <- c('cazyme' = '#46ACC8', 
            'lipase' = '#E2D200', 
            'peptidase' = '#E58601',
            'virulence' = "#B40F20")

## overall numbers
ggplot(all_secreted_summary, aes(x = genome, y = count, fill = type)) +
  geom_col() +
  labs(x = '', y = '') +
  theme_classic() +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave('plots/secretome_summary.pdf')





not_secreted <- all_enzymes %>% filter(!gene %in% secreted$gene)

not_secreted_summary <- not_secreted %>% 
  group_by(genome, type) %>%
  summarize(count = sum(count))

## overall numbers
ggplot(not_secreted_summary, aes(x = genome, y = count, fill = type)) +
  geom_col() +
  labs(x = '', y = '') +
  theme_classic() +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave('plots/non_secreted_catabolite_summary.pdf')



## just look at cazy num genes v num domains plot
ag_colours <- c("Generalist" = "#808080",
             "Coral" = "#CE3DD0",
             "Higher" = "#2D71F6",
             "Lower" = "#FFFEAB",
             "Leafcutter" = "#377D22",
             "Outgroup1" = 'black',
             "Outgroup2" = 'gray')

cazy_sum <- cazy %>% group_by(genome) %>%
  summarize(num_genes = length(unique(gene)),
            num_domain = length(unique(cazyme))) %>%
  rename('acc' = genome) %>%
  left_join(., metadata, by = 'acc') %>%
  mutate(Agriculture =ifelse(is.na(Agriculture), yes = 'Generalist', no = Agriculture))

cazy_sum %>%
  ggplot(., aes(x = num_genes, y = num_domain)) +
  geom_point(shape = 21, size = 2, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  scale_fill_manual(values = ag_colours) +
  theme_classic() +
  ggsave('plots/cazy_domain_diversity_all.pdf')




