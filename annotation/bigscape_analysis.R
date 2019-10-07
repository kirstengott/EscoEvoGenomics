source('~/scripts/theme_kirsten.R')
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
bigscape <- read_csv('bigscape.csv')






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

