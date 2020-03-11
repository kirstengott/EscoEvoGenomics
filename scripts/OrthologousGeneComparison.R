
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)
library(UpSetR)
library('ComplexHeatmap')

annot <- read_tsv('annotation/all_annotations.txt', col_names = c('acc', 'gene', 'tool', 'annot'))

meta <- read_csv('tables/metadata.csv')


all_gene_id_files <- list.files('annotation/fastortho/proteins/', pattern = 'gene_ids', full.names = TRUE)

all_gene_ids <- lapply(all_gene_id_files, function(x){
  ids <- scan(x, what = 'character')
  data.frame(gene = ids, genome = rep(basename(x), length(ids)), stringsAsFactors = FALSE) %>%
    mutate(acc = case_when(
      grepl('GCA_004303015', genome) ~ 'GCA_004303015.1',
      grepl('GCA', genome) ~ paste(strsplit(genome, split = "_")[[1]][c(1,2)], collapse = "_"),
      grepl("SPDT00000000.1_genomic", genome) ~ "SPDT00000000",
      grepl('LGSR', genome) ~ "LGSR00000000",
      grepl('.fasta_gene_ids', genome) ~ sub(".all.maker.proteins.fasta_gene_ids", '', genome),
      TRUE ~ genome
    )) %>%
    select(-genome) 
}) %>% bind_rows()

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
  full_join(all_gene_ids, ., by = c('gene', 'acc')) %>%
  left_join(., meta, by = 'acc') %>%
  mutate(orthology_meta = ifelse(is.na(OrthoGroup),
                                 yes = gene,
                                 no = OrthoGroup))


## compare the three pangenomes of the three clades



groups <- c('Trichoderma', 'Hypomyces/Cladobotryum', 'Escovopsis')

venn <- lapply(groups, function(x){
  if(x %in% c('Escovopsis', 'Trichoderma')){
    cutoff <- trunc(length(which(meta$clade_groups %in% x)) *.95)
  } else{
    cutoff <- length(which(meta$clade_groups %in% x))
  }
  ortho_full %>% 
    filter(clade_groups %in% x) %>%
    group_by(orthology_meta) %>%
    mutate(n_genomes = length(unique(genus_species))) %>%
    ungroup() %>%
    filter(n_genomes >= cutoff)
}) %>% bind_rows() %>%
  select(orthology_meta, clade_groups, n_genomes) %>%
  distinct() %>%
  group_by(clade_groups, orthology_meta) %>%
  summarize(genome_count = sum(n_genomes)) %>%
  ungroup()


esco_set <- venn %>% filter(clade_groups == 'Escovopsis') %>% .$orthology_meta %>% unique()
ch_set   <- venn %>% filter(clade_groups == 'Hypomyces/Cladobotryum') %>% .$orthology_meta %>% unique()
t_set    <- venn %>% filter(clade_groups == 'Trichoderma') %>% .$orthology_meta %>% unique()

listInput <- make_comb_mat(list('Escovopsis' = esco_set, 
                  'Hypomyces/Cladobotryum' = ch_set,
                  'Trichoderma' = t_set))
set_name(listInput)

pdf("plots/orthologues_upset.pdf", width = 8, height = 8)
 UpSet(t(listInput))#, order.by = "freq", 
#       point.size = 3.5, 
#       line.size = 1.5, 
#       text.scale = c(1.3, 1.3, 1, 1, 2, 0.75),
#       mb.ratio = c(0.8, 0.2),
#       empty.intersections = "on")
dev.off()
 
 ## Keep this as it show a combination of the set size of esco, with the intersection size with the others
data.frame(Escovopsis = length(esco_set),
           Overlap = length(intersect(which(esco_set %in% ch_set), which(esco_set %in% t_set)))) %>%
  gather(set, n_genes) %>%
  ggplot(aes(y = n_genes, x = set)) + geom_col(fill = 'black') +
  theme_minimal() + 
  labs(y = 'Number of Orthologous Gene Groups', 
       x = '') +
  ggsave('plots/esco_v_overlap_barchart.pdf', width = 5, height = 7)


venn.diagram(
  x = list('Escovopsis' = esco_set, 'Hypomyces/Cladobotryum' = ch_set, 'Trichoderma' = t_set),
  height = 3000,
  width = 3000,
   category.names = c(paste0("Escovopsis\n(", n_distinct(esco_set), " Orthologues across ", nrow(filter(meta, clade_groups == 'Escovopsis')), " genomes)"),
                      paste0("Cladobotryum_Hypomyces\n(", n_distinct(ch_set), " Orthogroups across ", nrow(filter(meta, clade_groups == 'Hypomyces/Cladobotryum')), " genomes)"), 
                      paste0("Trichoderma\n(", n_distinct(t_set), " Orthogroups across ", nrow(filter(meta, clade_groups == 'Trichoderma')), " genomes)")),
  filename = 'plots/shared_genes_venn_diagramm.png',
  resolution = 300,
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


