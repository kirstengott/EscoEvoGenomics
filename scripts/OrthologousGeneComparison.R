
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)

annot <- read_tsv('annotation/all_annotations.txt', col_names = c('acc', 'gene', 'tool', 'annot'))

meta <- read_csv('tables/metadata.csv')

ortho_full <- read_csv('annotation/fastortho/orthologues_full.csv', col_names = c('OrthoGroup', 'genome', 'gene')) %>%
  group_by(genome) %>%
  mutate(acc = case_when(
    grepl('GCA_004303015', genome) ~ 'GCA_004303015.1',
    grepl('GCA', genome) ~ paste(strsplit(genome, split = "_")[[1]][c(1,2)], collapse = "_"),
    grepl("SPDT00000000.1_genomic", genome) ~ "SPDT00000000",
    grepl('all.maker', genome) ~ sub(".all.maker.proteins", '', genome),
    grepl('LGSR', genome) ~ "LGSR00000000",
    TRUE ~ genome
  )) %>%
  ungroup() %>%
  select(-genome) %>%
  left_join(., meta, by = 'acc')


venn <- ortho_full %>% select(OrthoGroup, genus_species) %>%
  mutate(group = case_when(
    grepl('Trichoderma', genus_species) ~ 'Trichoderma',
    grepl('Clado', genus_species) ~ 'Cladobotryum_Hypomyces',
    grepl('Hypo', genus_species) ~ 'Cladobotryum_Hypomyces',
    TRUE ~ 'Escovopsis'
   )) #%>% 
  # group_by(OrthoGroup) %>%
  # mutate(number_groups = length(unique(group))) %>%
  # ungroup()

esco_set <- venn %>% filter(group == 'Escovopsis') %>% .$OrthoGroup
ch_set <- venn %>% filter(group == 'Cladobotryum_Hypomyces') %>% .$OrthoGroup
t_set <- venn %>% filter(group == 'Trichoderma') %>% .$OrthoGroup

sum(length(esco_set), length(ch_set), length(t_set))

myCol <- brewer.pal(3, "Set1")
venn.diagram(
  x = list(esco_set, ch_set, t_set),
  height = 3000,
  width = 3000,
  category.names = c(paste0("Escovopsis\n(", n_distinct(esco_set), " Orthologues across ", n_distinct(filter(venn, group == 'Escovopsis')%>% .$genus_species), " genomes)"),
                     paste0("Cladobotryum_Hypomyces\n(", n_distinct(ch_set), " Orthogroups across ", n_distinct(filter(venn, group == 'Cladobotryum_Hypomyces')%>% .$genus_species), " genomes)"), 
                     paste0("Trichoderma\n(", n_distinct(t_set), " Orthogroups across ", n_distinct(filter(venn, group == 'Trichoderma')%>% .$genus_species), " genomes)")),
  filename = 'plots/shared_genes_venn_diagramm.png',
  main = paste(n_distinct(venn$OrthoGroup), 'Orthogroups'),
  resolution = 300,
  # Circles
  fill = myCol,
  # Numbers
  fontface = "bold",
  fontfamily = "sans",
  main.fontface = 'bold',
  main.fontfamily = 'sans',
  # Set names
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)



ortho_annot <- ortho_full %>%
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


