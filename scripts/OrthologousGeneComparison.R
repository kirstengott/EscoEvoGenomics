
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

cazy <- read_csv('tables/cazy_annotation.csv') %>% 
  .$gene

card <- read_csv('tables/card_annotation.csv') %>%
  select(ORF_ID) %>%
  mutate(ORF_ID = sub(" .*$", "", ORF_ID)) %>%
  .$ORF_ID


length(which(ortho_full$gene %in% card))
length(which(ortho_full$gene %in% cazy))



ortho_full <- read_csv('annotation/fastortho/orthologues_full.csv', col_names = c('OrthoGroup', 'genome', 'gene')) %>%
  mutate(cazy = ifelse(gene %in% cazy,
                       yes = 'cazy',
                       no = NA)) %>%
  mutate(card = ifelse(gene %in% card,
                       yes = 'card',
                       no = NA)) %>%
  select(-gene, -genome) %>%
  filter(!is.na(cazy) | !is.na(card)) %>%
  group_by(OrthoGroup) %>%
  mutate(annotation = paste(na.omit(c(unique(cazy), unique(card))), collapse = "_")) %>%
  ungroup()  %>%
  select(-cazy, -card) %>%
  distinct() %>%
  data.frame(., stringsAsFactors = FALSE)
  

rownames(ortho_full) <- ortho_full$OrthoGroup
ortho_full$OrthoGroup <- NULL


ortho <- read_csv('annotation/fastortho/orthologues_presence_absence.csv')
ortho_mat <- data.frame(ortho, stringsAsFactors = FALSE)
rownames(ortho_mat) <- ortho$OrthoGroup
ortho_mat$OrthoGroup <- NULL
ortho_mat <- as.matrix(ortho_mat)


ortho_mat[ortho_mat > 1] <- 1


breaks = seq(0, 1)

#colors_f <- colorRampPalette(colors = c('white', 'black'))
colors <- c('white', 'black')


annot_colors <- list(annotation = c('card' =  "#1B9E77", 'cazy' = "#D95F02", "cazy_card" = "#7570B3"))

pheatmap(ortho_mat, show_rownames = FALSE, color = colors, cellwidth = 10, annotation_row = ortho_full, annotation_colors = annot_colors, 
         annotation_names_row = FALSE, filename = 'plots/orthologous_gene_matrix.pdf')


