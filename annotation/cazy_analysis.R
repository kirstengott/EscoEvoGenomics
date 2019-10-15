source('~/scripts/theme_kirsten.R')
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(ggfortify)

cazy_files <- list.files('cazy', full.names = TRUE)
# c <- cazy_files[1]
# 
# cazy <- lapply(cazy_files, function(c){
#   read_tsv(c) %>%
#     summarize(NumCazy = length(unique(`Gene ID`)))  %>%
#     mutate(genome = sub(".txt", "", basename(c)))
# }) %>% bind_rows() %>% as.data.frame()
# 
metadata <- read_csv('../metadata.csv') %>%
  mutate(acc_red = sub("\\..*$", "", acc))

# cazy$genome <- factor(cazy$genome, levels = rev(metadata$acc_red))
# ## cazy
# cazy %>% ggplot(., aes(y = NumCazy, x = genus_species)) +
#   geom_col() +
#   theme_kirsten(rot = 90) + coord_flip() +
#   ggsave('cazy.pdf')


#rownames(cazy) <- cazy$genome

#cazy[metadata$acc_red, 'NumCazy', drop = FALSE]



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

autoplot(pca_data, data = c_heat,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3) + ggsave('test.pdf')
