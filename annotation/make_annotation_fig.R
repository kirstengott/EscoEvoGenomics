source('~/scripts/theme_kirsten.R')
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
bigscape <- read_csv('bigscape.csv')


munge_gca <- function(x){
  paste(strsplit(x, split = "_")[[1]][c(1,2)], collapse = "_")
}

card_files <- list.files('card', full.names = TRUE)

card <- lapply(card_files, function(x){
  read_tsv(x) %>%
    rename('Resistance_Mechanism' = `Resistance Mechanism`) %>%
    select(ORF_ID, Resistance_Mechanism) %>%
    group_by(Resistance_Mechanism) %>%
    summarize(n_genes = length(unique(ORF_ID))) %>%
    mutate(genome = ifelse(grepl('GCA', x), yes = munge_gca(basename(x)),
                           no = sub("\\..*$", "", basename(x))))
}) %>% bind_rows()


cazy_files <- list.files('cazy', full.names = TRUE)
c <- cazy_files[1]

cazy <- lapply(cazy_files, function(c){
  read_tsv(c) %>%
    summarize(NumCazy = length(unique(`Gene ID`)))  %>%
    mutate(genome = sub(".txt", "", basename(c)))
}) %>% bind_rows() %>% as.data.frame()

metadata <- read_csv('../Figures/bigscape/metadata.csv') %>%
  mutate(acc_red = sub("\\..*$", "", acc))

cazy$genome <- factor(cazy$genome, levels = rev(metadata$acc_red))
## cazy
cazy %>% ggplot(., aes(y = NumCazy, x = genome)) +
  geom_col() +
  theme_kirsten(rot = 90) + coord_flip() +
  ggsave('cazy.pdf')


rownames(cazy) <- cazy$genome

cazy[metadata$acc_red, 'NumCazy', drop = FALSE]

## bigscape heatmap
big_mat <- bigscape %>% spread(Type, total) %>% as.data.frame()
rownames(big_mat) <- big_mat$Genome
big_mat$Genome <- NULL
big_mat <- as.matrix(big_mat)
colors <- brewer.pal(4, 'Purples')
pheatmap(big_mat[metadata$acc,],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 20,
         cellheight = 20,
         color = colors, filename = 'big_heat.pdf')



card_mat <- card %>% filter(n_genes >= 10) %>%

  spread(Resistance_Mechanism, n_genes) %>% as.data.frame()
rownames(card_mat) <- card_mat$genome
card_mat$genome <- NULL
card_mat <- as.matrix(card_mat)
#colors <- brewer.pal(11, 'Blues')
colors <- colorRampPalette(brewer.pal(n = 7, name ="Blues"))(11)
card %>% filter(n_genes >= 10) %>% .$n_genes %>% summary()
breaks <- seq(min(card$n_genes), max(card$n_genes), by = 30)


pheatmap(card_mat[metadata$acc,],
                    breaks = breaks,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     cellwidth = 20,
                     cellheight = 20,
                     color = colors, filename = 'card_heat.pdf')

head(card)
head(bigscape)
head(cazy)

fam_db <- read_tsv('CAZyDB.07312018.fam-activities.txt', comment = '#', col_names = c('cazy_base', 'cazy_description'))

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
  mutate(n_genes = as.numeric(length(unique(gene))))

c_sum <- cazy %>% ungroup() %>% select(-start, -stop, -gene) %>% distinct() %>% spread(genome, n_genes)

c_heat <- c_sum %>% select(-cazy_base, -cazy_description) %>% data.frame()
c_heat[is.na(c_heat)] <- 0

rownames(c_heat) <- c_heat$cazy
c_heat$cazy <- NULL


c_heat_f <- c_heat[which(rowSums(c_heat) > 20), ]

colors <- colorRampPalette(brewer.pal(n = 7, name ="Reds"))(11)
pheatmap(as.matrix(c_heat),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border_color = "grey70",
         cellwidth = 10,
         cellheight = 10, 
         filename = 'cazy_heat.pdf',
         color = colors)

pca_data <- prcomp(c_heat)

metadata_t <- metadata %>% mutate(acc = sub("\\..*$", "", acc))

pca_rot <- data.frame(pca_data$rotation,
                      acc = rownames(pca_data$rotation)
) %>%
  left_join(., metadata_t, by = 'acc') %>%
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

ggplot(data = pca_rot, aes(x = PC1, y = PC2)) +
  geom_point(shape = 21, size = 4, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  #geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 5)
  theme_classic() +
  #geom_text(data = geom_text_data, hjust = 1, vjust = 0, size = 3) +
  scale_fill_manual(values = colours) +
  scale_x_continuous(expand = c(0, 0.2)) +
  ggsave('cazy_pca_12.pdf')

eigs <- pca_data$sdev^2
proportion = (eigs/sum(eigs))*100
cumulative = cumsum(eigs)/sum(eigs)
## look at this
head(pca_data$x)
## bigscape PCA
library(ggfortify)

autoplot(pca_data, data = c_heat,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3) + ggsave('test.pdf')
 

shorten_genus_species <- function(x){
  spl = strsplit(x, split = " ")
  paste0(strsplit(spl[[1]][1], split = "")[[1]][1], '. ', spl[[1]][2])
}

metadata <- read_csv('../Figures/bigscape/metadata.csv') %>%
  mutate(acc_red = sub("\\..*$", "", acc)) %>%
  group_by(genus_species) %>%
  mutate(pca_label = ifelse(is.na(Agriculture),
                            yes = shorten_genus_species(genus_species),
                            no = genus_species)) %>%
  ungroup()


b_pca <- bigscape %>% spread(key = Genome, value = total) %>% data.frame(.)
rownames(b_pca) <- b_pca$Type
b_pca$Type <- NULL

b_pca[is.na(b_pca)] <- 0

pca_data <- prcomp(b_pca)


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
  geom_point(shape = 21, size = 4, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  #geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 5)
  theme_classic() +
  geom_text(data = geom_text_data, hjust = 1, vjust = 0, size = 3) +
  scale_fill_manual(values = colours) +
  scale_x_continuous(expand = c(0, 0.2)) +
  ggsave('bigscape_pca_12.pdf')

ggplot(data = pca_rot, aes(x = PC1, y = PC3, label = pca_label)) +
  geom_point(shape = 21, size = 4, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  theme_classic() +
  geom_text(data = geom_text_data, hjust = 1, vjust = 0, size = 3) +
  scale_colour_manual(values = colours) +
  scale_x_continuous(expand = c(0, 0.2)) +
  ggsave('bigscape_pca_13.pdf')



## PCA without trichodermat
keep_accessions <- metadata %>% filter(genus != "Trichoderma") %>% .$acc %>% unique()

b_pca <- bigscape %>% filter(Genome %in% keep_accessions) %>%
  spread(key = Genome, value = total) %>% data.frame(.)
rownames(b_pca) <- b_pca$Type
b_pca$Type <- NULL

b_pca[is.na(b_pca)] <- 0

pca_data <- prcomp(b_pca)


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
  theme_classic() +
  geom_text(data = geom_text_data, hjust = 1, vjust = 0, size = 3) +
  scale_fill_manual(values = colours) +
  scale_x_continuous(expand = c(0, 0.2)) +
  ggsave('bigscape_pca_12_notrich.pdf')


ggplot(data = pca_rot, aes(x = PC1, y = PC3, label = pca_label)) +
  geom_point(shape = 21, size = 4, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  theme_classic() +
  geom_text(data = geom_text_data, hjust = 1, vjust = 0, size = 3) +
  scale_fill_manual(values = colours) +
  scale_x_continuous(expand = c(0, 0.2)) +
  ggsave('bigscape_pca_13_notrich.pdf')


ggplot(data = pca_rot, aes(x = PC2, y = PC3, label = pca_label)) +
  geom_point(shape = 21, size = 4, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  theme_classic() +
  geom_text(data = geom_text_data, hjust = 1, vjust = 0, size = 3) +
  scale_fill_manual(values = colours) +
  scale_x_continuous(expand = c(0, 0.2)) +
  ggsave('bigscape_pca_23_notrich.pdf')

## for CAFE
# The data file may consist of data on one family or up to thousands of families for the specified tree.
# The first line of the data file should contain the extant species’ names (as used in the Newick tree description),
# tab-delimited in no particular order.
# Subsequent lines each correspond to a gene family and contain tab-delimited family sizes for these extant species.
# Columns in the data file whose header does not correspond to any of the names in the Newick tree description,
# which may provide additional information about the gene families, are ignored and simply copied in the output file


cazy_caf <- cazy %>% select(genome, cazy_base, cazy_description, n_genes) %>% distinct() %>%
#  spread(key = genome, value = n_genes) %>%
  rename('family' = cazy_base,
         'family_description' = cazy_description)
glimpse(cazy_caf)

bigscape_caf <- bigscape %>% mutate(Genome = sub("\\..*$", "", Genome)) %>%
  #spread(key = Genome, value = total) %>%
  rename('family' = Type,
         'genome' = Genome,
         'n_genes' = total) %>%
  mutate(family_description = 'null') #%>% select(family, family_description, everything())

glimpse(bigscape_caf)

card_caf <- card %>% mutate(genome = sub("\\..*$", "", genome)) %>%
  #spread(key = genome, value = n_genes) %>%
  rename('family' = Resistance_Mechanism) %>%
  mutate(family_description = 'null') #%>% select(family, family_description, everything())

glimpse(card_caf)

metadata <- read_tsv('id_metadata.txt') %>%
  mutate(genome = sub("\\..*$", "", acc)) %>%
  select(-acc, -family, -genus, -genome_id, -id)

all_cafe <- bind_rows(cazy_caf, bigscape_caf, card_caf) %>%
  ungroup() %>%
  left_join(., metadata, by = 'genome') %>%
  select(-genome) %>%
  spread(key = genus_species, value = n_genes)


all_cafe[is.na(all_cafe)] <- 0


glimpse(all_cafe)

write_tsv(all_cafe, path = 'cafe_table_input.tsv')


# ranked_cazy <-
#   cazy %>% select(cazy_base, cazy_description, n_genes) %>% distinct() %>% ungroup() %>% mutate(rank = rank(-n_genes))
#
# ranked_cazy %>% filter(rank <= 10)
# cazy

