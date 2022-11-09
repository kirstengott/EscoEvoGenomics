library(tidyverse)
library(cowplot)

## GCA_008477525.1 = SPDT00000000
source('scripts/working_genomes.R') 

genomes_keep <- genomes_keep %>% sub(., pattern = '\\.1', replacement = '')

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
           'GCA_008477525',
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

order = order[order %in% genomes_keep]

meta <- read_csv('tables/metadata.csv')  %>% select(acc, Agriculture) %>% na.omit() %>%
  mutate(g1 = sub("\\..*$", "", acc)) %>%
  select(-acc) %>%
  rename("Agriculture_g1" = Agriculture)

meta2 <- meta %>% 
  mutate(g2 = g1) %>% 
  select(-g1) %>%
  rename("Agriculture_g2" = Agriculture_g1)


dn_ds <- read_csv('data/mafft-dn_ds_all.txt', 
                  col_names = c('orthogroup', 'g1', 'g2', 'ds', 'dn', 'dn_ds')) %>%
  mutate(g1 = sub("\\..*$", "", g1)) %>%
  mutate(g2 = sub("\\..*$", "", g2)) %>%
  filter(g1 %in% order, g2 %in% order) %>%  
  left_join(., meta, by = 'g1') %>%
  left_join(., meta2, by = 'g2')

dn_ds %>%
  filter(dn < 10, ds < 10) %>% summary()

dn_ds %>%
  filter(dn < 10, ds < 10, dn_ds < 10) %>%
  ggplot(aes(y = dn_ds, x = g1)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dn_ds %>%  
  filter(ds <= 10) %>%
  ggplot(aes(x = ds)) + geom_histogram()

dn_ds %>%
  filter(dn <= 10) %>%
  ggplot(aes(x = dn)) + geom_histogram()

source('scripts/color_palettes.R')



## only examining dn/ds within lineages 

df <- dn_ds %>% filter(Agriculture_g1 == Agriculture_g2) %>%
  mutate(Agriculture = Agriculture_g1) %>%
  filter(dn < 10, ds < 10, dn_ds < 5)

df$Agriculture <- factor(df$Agriculture, levels = c("Trichoderma", 
                                                    "Hypomyces_Cladobotryum",
                                                    "Lower", 
                                                    "Coral", 
                                                    "Higher", "Leafcutter"))
                               


  
a <- df %>%  
  ggplot(., aes(x = Agriculture, y = dn_ds, fill = Agriculture)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = colors) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+
#  labs(title = ,
#       subtitle = )


paste0('Pairwise dN/dS for ', n_distinct(df$orthogroup), " Orthologues")
paste0(nrow(df), ' total observations')


load('rdata/orthogroup_pa_ag.rdata')

lin_spec_genes <- orthogroup_pa_ag %>% 
  filter(lineage_specific == TRUE)

lin_spec <- df %>% group_by(orthogroup) %>%
  mutate(lineage_specific = ifelse(orthogroup %in% lin_spec_genes$orthogroup, 
                                   yes = "Lineage Specific", 
                                   no = "Shared Across Lineages"))



nrow(lin_spec %>% filter(lineage_specific == "Lineage Specific"))
nrow(lin_spec %>% filter(lineage_specific == "Shared Across Lineages"))

lin_spec$lineage_specific <- factor(lin_spec$lineage_specific, 
                                    levels = c("Shared Across Lineages", "Lineage Specific"))

b <- lin_spec %>%
  ggplot(., aes(y = dn_ds, x = Agriculture, fill = Agriculture)) +
  geom_boxplot() +
  scale_fill_manual(values = colors) + 
  facet_grid(~lineage_specific) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray65')

b + ggsave(filename = 'plots/dn_ds.pdf')


# c <- plot_grid(plotlist = list(a, b),
#           nrow = 2, labels = 'AUTO')
#   
# cowplot::save_plot(c, filename = 'plots/dn_ds.pdf', base_height = 7,
#                   base_width = 5)
  







genomes_keep <- genomes_keep %>% sub(., pattern = '\\.1', replacement = '')

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
           'GCA_008477525',
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

order = order[order %in% genomes_keep]




data <- read_tsv('tables/ANIvis.tsv', col_types = 'ccnn') %>% 
  mutate_at(c('Query', 'Reference'), basename) %>%
  mutate_at(c('Query', 'Reference'), function(x){sub(".cds_from_genomic.fa", "", x)}) %>%
  mutate_at(c('Query', 'Reference'), function(x){sub("\\..*$", "", x)}) %>%
#  select(-QueryAligned) %>% 
  rename('acc' = Query) %>% 
  mutate(acc = ifelse(acc == "SPDT00000000", yes = "GCA_008477525", no = acc)) %>%
  mutate(Reference = ifelse(Reference == "SPDT00000000", yes = "GCA_008477525", no = Reference)) %>%
  filter(acc %in% genomes_keep, Reference %in% genomes_keep) %>%
  rename('g1' = acc, 'g2' = Reference) %>%
  rename('ani' = ANI)

data %>% 
  rename('ANI' = ani) %>%
  ggplot(aes(x = ANI, y = QueryAligned)) +
  geom_point() +
  theme_classic() +
  geom_vline(xintercept = 98, linetype = "dashed", color = 'gray37') +
  theme(legend.position = 'none') +
  ggsave('plots/ani_queryaligned.pdf')


data <- data %>% select(-QueryAligned)








#t <- order[10]
datalist = list()
data_loop <- data
for (t in order) {
  print(t)
  ani_vals <- c(data_loop[grep(t, data_loop$g1), ]$ani,
                data_loop[grep(t, data_loop$g2), ]$ani)
  
  genome_vals <- c(data_loop[grep(t, data_loop$g1), 'g2']$g2, 
                   data_loop[grep(t, data_loop$g2), 'g1']$g1)
  
  #data_t <- data_loop %>% filter(g1 != t, g2 != t)
  df_out  <- data.frame(g1 = t, 
                        g2 = genome_vals, 
                        ani = ani_vals) %>%
    rbind(., data.frame(g1 = t, g2 = t, ani = 100))
  data_loop <- data_loop %>% filter(g2 != t)
  df_out1 <- lapply(order, function(x){ ## reorder the df
    filter(df_out, g2 == x)
  }) %>% bind_rows() %>% group_by(g1, g2) %>%
    mutate(ani = mean(ani)) %>% distinct()
  
  datalist[[t]] <- df_out1
}


ani <- bind_rows(datalist)
ani$g1 <- factor(ani$g1, levels = order)
ani$g2 <- factor(ani$g2, levels = order)


library(RColorBrewer)
palette <- c(colorRampPalette(colors = c('white', 'gray37'))(95), 'black')




ggplot(data = ani, aes(y = g1, x = g2, fill = ani)) +
  geom_tile() + 
  theme_classic() +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = 'left') +
  theme(axis.text.x = element_text(angle = 65, hjust = 0)) +
  scale_fill_steps(breaks = c(seq(1, 98), 100), low = 'white', high = 'black') +
  labs(y = '', x = '') +
  ggsave('plots/ani.pdf')




data$g1 <- factor(data$g1, levels = order)
data$g2 <- factor(data$g2, levels = order)


ggplot(data = data, aes(y = g1, x = g2, fill = ani)) +
  geom_tile() + 
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = 'left') +
  theme(axis.text.x = element_text(angle = 65, hjust = 0)) +
  scale_fill_steps(breaks = c(seq(1, 98), 100), low = 'white', high = 'black') +
  labs(y = '', x = '') #+
#  ggsave('plots/ani.pdf')




ani_groups <- c(rep(1, 8), 2, 2, 3, 3, 4, 4, 5, 5, rep(7, 3),
                8, 8, 8, 9, 9, 9, 9, 10, 10, 11, 11, 11)


ani_groups <- cbind(order, ani_groups) %>% data.frame() %>% select(order, ani_groups)

write_csv(ani_groups, path = 'tables/ani_groups.csv')








