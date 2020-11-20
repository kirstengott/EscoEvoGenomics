library(tidyverse)


data <- read_csv('tables/ani_dn_ds_summary.txt', 
                 col_names = c(
  'g1', 'g2', 'dn_ds', 'ani'
))


order <- c("GCA_000225605",
           "GCA_003012105",
           "GCA_003025115",
           "GCA_000171015",
           "GCA_001481775",
           "GCA_002894205",
           "GCA_003025105",
           "GCA_003025155",
           "GCA_001050175",
           "GCA_000167675",
           "GCA_000513815",
           "GCA_000170995",
           "GCA_002894145",
           "GCA_000988865",
           "GCA_011066345",
           "GCA_002022785",
           "GCA_003025095",
           "GCA_002838845",
           "SPDT00000000",
           "GCA_004303015",
           "GCA_011799845",
           "ICBG2046",
           "ICBG2048",
           "ICBG2047",
           "ICBG2049",
           "ICBG712",
           "ICBG721",
           "ICBG1054",
           "ICBG726",
           "ICBG1065",
           "ICBG1075",
           "NIGD00000000",
           "ICBG710",
           "ICBG730",
           "ICBG751",
           "ICBG733",
           "ICBG1096",
           "ICBG742",
           "NIGB00000000",
           "LGSR00000000",
           "NQYS00000000",
           "ICBG731",
           "ICBG736",
           "NIGC00000000",
           "NQYR00000000",
           "NQYQ00000000")


datalist = list()
for (t in order) {
  print(t)
  ani_vals <- c(data[grep(t, data$g1), ]$ani,
                data[grep(t, data$g2), ]$ani)
  genome_vals <- c(data[grep(t, data$g1), ]$g2, 
                   data[grep(t, data$g2), ]$g1)
  data <- data %>% filter(g1 != t, g2 != t)
  datalist[[t]] <- data.frame(g1 = t, 
                              g2 = genome_vals, 
                              ani = ani_vals)
}


ani <- bind_rows(datalist)







data <- read_csv('tables/ani_dn_ds_summary.txt', 
                 col_names = c(
                   'g1', 'g2', 'dn_ds', 'ani'
                 ))

datalist = list()
for (t in order) {
  print(t)
  ani_vals <- c(data[grep(t, data$g1), ]$dn_ds,
                data[grep(t, data$g2), ]$dn_ds)
  genome_vals <- c(data[grep(t, data$g1), ]$g2, 
                   data[grep(t, data$g2), ]$g1)
  data <- data %>% filter(g1 != t, g2 != t)
  datalist[[t]] <- data.frame(g1 = t, 
                              g2 = genome_vals, 
                              dn_ds = ani_vals)
}


dn_ds <- bind_rows(datalist)

data <- left_join(dn_ds, ani, by = c('g1', 'g2'))

data$g1 <- factor(data$g1, levels = order)
data$g2 <- factor(data$g2, levels = order)

data %>% 
  ggplot() +
  geom_tile(aes(y = g2, x = g1, fill = ani)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggsave(filename = 'plots/ani_plot.pdf', height = 6,
         width = 7)

data %>% 
  ggplot() + 
  geom_tile(aes(y = g1, x = g2, fill = dn_ds)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_continuous(type = 'viridis')
  ggsave(filename = 'plots/dn_ds_plot.pdf', height = 6,
         width = 7)


data <- read_csv('tables/mafft-dn_ds_summary.txt', col_names = c('g1', 'g2', 'dn_ds')) %>%
  mutate(g1 = sub("\\..*$", "", g1)) %>%
  mutate(g2 = sub("\\..*$", "", g2))
 

datalist = list()
for (t in order) {
  print(t)
  ani_vals <- c(data[grep(t, data$g1), ]$dn_ds,
                data[grep(t, data$g2), ]$dn_ds)
  genome_vals <- c(data[grep(t, data$g1), ]$g2, 
                   data[grep(t, data$g2), ]$g1)
  data <- data %>% filter(g1 != t, g2 != t)
  datalist[[t]] <- data.frame(g1 = t, 
                              g2 = genome_vals, 
                              dn_ds = ani_vals)
}


data <- bind_rows(datalist)

#data <- left_join(dn_ds, ani, by = c('g1', 'g2'))

data$g1 <- factor(data$g1, levels = order)
data$g2 <- factor(data$g2, levels = order)





data %>% 
  ggplot() + 
  geom_tile(aes(y = g1, x = g2, fill = dn_ds)) +
  theme_classic() +
  scale_fill_continuous( type = 'viridis') +
  theme(axis.text.x = element_text(angle = 90))

##################################
### PLAYING WITH DENSITY PLOTS ###
##################################



metadata <- read_csv('tables/metadata.csv') %>%
  select(acc, genus, Agriculture) %>%
   filter(!is.na(genus)) #%>%
  # mutate(acc = sub("\\..*$", "", acc)) %>%
  # mutate(Query = acc) %>%
  # mutate(Reference = acc)

esc <- metadata %>% filter(genus == 'Escovopsis') %>% .$acc
trich <- metadata %>% filter(genus == 'Trichoderma') %>% .$acc %>%
  purrr::map(sub, pattern = "\\..*$", replacement = "") %>%
  unlist()

data %>% 
  filter(g1 %in% esc, g2 %in% esc) %>%
  ggplot(aes(x = dn_ds, y = ani)) +
  geom_density_2d() +
  geom_point() +
  theme_classic()

data %>% 
    filter(g1 %in% trich, g2 %in% trich) %>%
  ggplot(aes(x = dn_ds, y = ani)) +
    geom_density_2d() +
    geom_point() +
    theme_classic()



split2vec <- function(x, split){
  paste(strsplit(x, split = split)[[1]][c(1,2)], collapse = "_")
}

ani <- read_tsv('tables/ANIvis.tsv') %>%
  mutate_at(c('Query', 'Reference'), basename) %>%
  mutate_at(c('Query', 'Reference'), sub, pattern = '.cds_from_genomic.fa', replacement = '')

ani$Query <- ani$Query %>% purrr::map(split2vec, split = '_') %>% 
  purrr::map(sub, pattern = '_NA', replacement ='') %>% 
  unlist

ani$Reference <- ani$Reference %>% purrr::map(split2vec, split = '_') %>% 
  purrr::map(sub, pattern = '_NA', replacement ='') %>% 
  unlist



ani_d <- ani %>% 
  left_join(., select(metadata, -Reference), by = 'Query') %>%
  left_join(., select(metadata, -Query), by = 'Reference') %>%
  mutate(genus_join = paste0(genus.x, '-', genus.y)) %>%
  mutate(group_join = paste0(Agriculture.x, '-', Agriculture.y))

poin_sub <- ani %>%
  filter(ANI > 95)

ani %>% 
  filter(Query %in% esc, Reference %in% esc) %>%
  ggplot(aes(x = QueryAligned, y = ANI)) +
  geom_point() +
  stat_density_2d() + 
  theme_classic()

ani %>% 
  filter(Query %in% trich, Reference %in% trich) %>%
  ggplot(aes(x = QueryAligned, y = ANI)) +
  geom_point() +
  stat_density_2d() + 
  theme_classic()

