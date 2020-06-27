library(pheatmap)
library(RColorBrewer)
library(pvclust)
library(vegan)
library(UpSetR)
source('scripts/color_palettes.R')
library(cowplot)
library(ggridges)
library(tidyverse)


make_ordered_dataframe <- function(data, order, col){
  datalist = list()
  for (t in order) {
    ani_vals <- c(data[grep(t, data$g1), col, drop = TRUE],
                  data[grep(t, data$g2), col, drop = TRUE])
    genome_vals <- c(data[grep(t, data$g1), ]$g2, 
                     data[grep(t, data$g2), ]$g1)
    if(length(ani_vals > 0)){
      data <- data %>% filter(g1 != t, g2 != t)
      datalist[[t]] <- data.frame(g1 = t, 
                                  g2 = genome_vals, 
                                  dn_ds = ani_vals)
    } else {
      datalist[[t]] <- data.frame()
    }
  }
  data <- bind_rows(datalist)
  data$g1 <- factor(data$g1, levels = order)
  data$g2 <- factor(data$g2, levels = order)
  invisible(data)  
}

# dn_ds_f     <- read_csv('tables/mafft-dn_ds_all.txt', col_names = c('Orthogroup', 'g1', 'g2', 'dn_ds', 'dn', 'ds')) %>%
#   mutate(g1 = sub("\\..*$", "", g1)) %>%
#   mutate(g2 = sub("\\..*$", "", g2)) %>%
#   filter(dn != 99, 
#          ds != 99)

esc <- c(           "ICBG2046",
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


#dn_ds_f$g1 <- factor(dn_ds_f$g1, levels = order)
#dn_ds_f$g2 <- factor(dn_ds_f$g2, levels = order)

load('rdata/orthologues_long.rdata')



genomes <- colnames(orthog) %>% 
  purrr::map(.x = ., function(x){sub(pattern = ".all.maker.proteins", replacement = "", x)}) %>%
  purrr::map(.x = ., function(x){sub(pattern = ".protein.", replacement = "", x)}) %>%
  purrr::map(.x = ., function(x){sub(pattern = ".protein", replacement = "", x)}) %>%
  purrr::map(.x = ., function(x){split <- strsplit(x, split = "_")[[1]]
                                                   if (length(split) > 1) {
                                                     paste(split[c(1,2)], collapse = "_")
                                                   } else {x} 
                                          
               }) %>%
  unlist()
genomes <- genomes[-c(1)]


# busco_f <-list.files('annotation/ascomybocta_odb10_busco/prot', full.names = TRUE)
# index <- sapply(genomes, function(x){
#   grep(x, busco_f)
# })

# 
# busco <- lapply(seq(1, length(index)), function(x){
#    genome = names(index[x])
#    df = data.frame()
#    try(
#    df <- read_tsv(busco_f[index[[x]]], skip = 3, col_names = c('busco', 'status',	'sequence', 'score', 	'length')) %>%
#      mutate(genome = genome))
#   return(df)
# }) %>% bind_rows() %>%
#    filter(!is.na(sequence))

# card_f <- list.files('annotation/card', full.names = TRUE)
# 
# index <- sapply(genomes, function(x){
#   grep(x, card_f)
# 
# #x <- 4
# card <- lapply(seq(1, length(index)), function(x){
#   genome = names(index[x])
#   df <- data.frame()
#   try(
#   df <- read_tsv(card_f[index[[x]]]) %>%
#     mutate(genome = genome))
#   return(df)
# }) %>% bind_rows() %>%
#   mutate(ORF_ID = sub(" .*$", "", ORF_ID))


metadata <- read_csv('tables/metadata.csv') %>% mutate(genome = sub("\\..*$", "", acc)) %>%
  mutate(Agriculture = ifelse(!clade_groups == 'Escovopsis', yes = clade_groups, no = Agriculture)) %>%
  mutate(Agriculture1 = ifelse(!clade_groups == 'Escovopsis', yes = clade_groups, no = Agriculture1))



orthog_meta <- orthog_long %>%
  group_by(Orthogroup) %>%
  mutate(n_genes = length(unique(genes))) %>%
  ungroup() %>%
  group_by(genome, Orthogroup) %>%
  mutate(n_genes_genome = length(unique(genes))) %>%
  filter(n_genes_genome <= 10) %>%
  ungroup() %>%
  mutate(g1 = sub("\\..*$", "", genome)) %>%
  mutate(genome = sub("\\..*$", '', genome)) %>%
  left_join(., metadata, by = 'genome') %>%
  filter(genome != 'GCA_000225605')

orthog_meta %>% group_by(clade_groups) %>%
  summarize(mean_copy_number = mean(n_genes_genome))

orthog_meta %>% 
    select(Orthogroup, genome, clade_groups, n_genes_genome) %>%
    distinct() %>%
    ggplot(aes(x = clade_groups, y = n_genes_genome)) +
    geom_violin() +
    theme_classic() +
    labs(y = 'Gene Copy Number', 
         title = 'Copy Number Distribution')
  

tool_sub <- c('virulence', 'led', 'MEROPS', 'cazy', 'CARD', 'busco', 'antismash')
tools_name <- c('virulence' = 'Virulence Genes', 
                'led' = 'Lipase Genes', 
                'MEROPS' = 'Peptidase Genes', 
                'cazy' = 'CAZyme Genes', 
                'CARD' = 'Resistance Genes', 
                'busco' = 'Single Copy Genes',
                'antismash' = 'BGC')


merops_fams <- read_tsv('annotation/merops_family_map.txt', col_names = c('annotation', 'merops_family')) %>%
  distinct()

merops_fam_ls        <- merops_fams$merops_family
names(merops_fam_ls) <- merops_fams$annotation

annot    <- read_tsv('annotation/all_annotations.txt', col_names = c('genome', 'gene', 'tool', 'annotation')) %>%
  filter(tool %in% tool_sub, annotation != '-') %>%
  mutate(genome = sub("\\..*$", "", genome)) %>%
  mutate(annotation = ifelse(tool == 'MEROPS',
                             yes = sub("-.*$", "", annotation),
                             no = annotation)) %>%
  mutate(annotation = ifelse(tool == 'MEROPS',
                             yes = ifelse(annotation %in% names(merops_fam_ls), yes = merops_fam_ls[annotation], no = annotation), ## integrate family ids, but keep undefined families
                             no = annotation)) %>%
  filter(!grepl('mix.component', annotation)) %>%
  mutate(annotation = ifelse(tool == 'virulence',
                             yes = sub("_.*$", "", annotation),
                             no = annotation))



orthog_long <- orthog_long %>%
  mutate(genome = sub("\\..*$", "", genome)) %>%
  filter(genome != "GCA_000225605") %>% ## remove cordycep genome for this analysis 
  rename('gene' = genes)


orthog_meta_ls <- lapply(tool_sub, function(x){
  if(x == 'led'){
    annot_sub <- annot %>% filter(tool == x) %>%
      group_by(genome, gene) %>%
      mutate(annotation = list(strsplit(annotation, split = ',')[[1]])) %>%
      ungroup() %>%
      unnest_longer(annotation)
  } else {
    annot_sub <- annot %>% filter(tool == x)}
  ## 90% or higher of the genes in the orthogroup need to be annotated to one of the tools to be included 
  orthog_meta <- orthog_long %>%
    inner_join(., annot_sub, by = c('genome', 'gene')) %>%
    mutate(metadata = ifelse(gene %in% annot_sub$gene, yes = TRUE, no = FALSE)) %>%
    group_by(Orthogroup) %>%
    mutate(n_genes = length(unique(gene))) %>%
    mutate(n_genes_in = length(which(metadata))) %>%
    ungroup() %>%
    group_by(genome, Orthogroup) %>%
    mutate(n_genes_genome = length(unique(gene))) %>%
    ungroup() %>%
    mutate(prop_genes_in = (n_genes_in/n_genes)*100) %>%
    mutate(g1 = sub("\\..*$", "", genome)) %>%
    filter(prop_genes_in >= 90) %>%
    mutate(genome = sub("\\..*$", '', genome)) %>%
    left_join(., metadata, by = 'genome') %>%
    mutate(tool = tools_name[x])
}) %>% bind_rows()


orthog_meta_ls %>% group_by(clade_groups, tool) %>%
    summarize(mean_copy_number = mean(n_genes_genome)) %>%
  spread(tool, mean_copy_number)

orthog_meta_ls$tool <- factor(orthog_meta_ls$tool, levels = c("Single Copy Genes",
                                                              "Resistance Genes",
                                                              "Virulence Genes",
                                                              "BGC",
                                                              "Lipase Genes",
                                                              "Peptidase Genes",
                                                              "CAZyme Genes"))

####################################################################################
## COPY NUMBER is the main reduction 
####################################################################################
#x <- tool_sub[1]

orthog_meta_ls %>% 
    select(Orthogroup, genome, clade_groups, n_genes_genome, tool) %>%
    distinct() %>%
    ggplot(aes(y = clade_groups, x = n_genes_genome, fill = clade_groups, color = clade_groups)) +
    geom_density_ridges(scale = .9, stat = 'binline', draw_baseline = FALSE, binwidth = 1) +
     scale_fill_manual(values = clade_colors) +
     scale_color_manual(values = clade_colors) +
    scale_y_discrete(expand = expansion(mult = c(0.01, .7))) +
    theme_classic() +
    theme(#is.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = 'none') +
    labs(x = 'Gene Copy Number',
         y = '') +
  facet_wrap(~tool, nrow = 3, ncol = 3) +
  ggsave('plots/gene_copy_number_clade_hist.pdf')


orthog_meta_ls %>% 
  select(Orthogroup, genome, clade_groups, n_genes_genome, tool) %>%
  distinct() %>%
  ggplot(aes(x = clade_groups, y = n_genes_genome, fill = clade_groups, color = clade_groups)) +
  geom_violin() +
  scale_fill_manual(values = clade_colors) +
  scale_color_manual(values = clade_colors) +
  # scale_y_discrete(expand = expansion(mult = c(0.01, .7))) +
  theme_classic() +
  theme(#axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = 'none') +
  labs(y = 'Gene Copy Number',
       x = '') +
  facet_wrap(~tool, nrow = 3, ncol = 3) +
  ggsave('plots/gene_copy_number_clade_violin.pdf')



orthog_meta_ls$Agriculture <- factor(orthog_meta_ls$Agriculture, levels = rev(c("Trichoderma",
                                                                            "Hypomyces_Cladobotryum",
                                                                            "Lower",
                                                                            "Coral",
                                                                            "Higher",
                                                                            "Leafcutter"
                                                                            )))

orthog_meta_ls %>% 
  select(Orthogroup, genome, Agriculture, n_genes_genome, tool) %>%
  distinct() %>%
  ggplot(aes(y = Agriculture, x = n_genes_genome, fill = Agriculture, color = Agriculture)) +
  geom_density_ridges(scale = .9, stat = 'binline', draw_baseline = FALSE, binwidth = 1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_y_discrete(expand = expansion(mult = c(0.01, .2))) +
  scale_x_continuous(breaks = seq(1, 30, by = 5), 
                     labels = seq(1, 30, by = 5),
                     expand = c(0, 0)) +
  theme_classic() +
  theme(#axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = 'none') +
  labs(x = 'Gene Copy Number',
       y = '') +
  facet_wrap(~tool, nrow = 3, ncol = 3) +
  ggsave('plots/gene_copy_number_ag_hist.pdf')


orthog_meta_ls$Agriculture <- factor(orthog_meta_ls$Agriculture, levels = c("Trichoderma",
                                                                                "Hypomyces_Cladobotryum",
                                                                                "Lower",
                                                                                "Coral",
                                                                                "Higher",
                                                                                "Leafcutter"
))

orthog_meta_ls %>% 
  select(Orthogroup, genome, Agriculture, n_genes_genome, tool) %>%
  distinct() %>%
  ggplot(aes(x = Agriculture, y = n_genes_genome, fill = Agriculture, color = Agriculture)) +
  geom_violin() +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  # scale_y_discrete(expand = expansion(mult = c(0.01, .7))) +
  theme_classic() +
  theme(#axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = 'none') +
  labs(x = '',
       y = 'Gene Copy Number') +
  facet_wrap(~tool, nrow = 3, ncol = 3) + 
  ggsave('plots/gene_copy_number_ag_violin.pdf')





###########################
## Making scatter plots ###
###########################



# df <- orthog_meta_ls %>%
#   select(Orthogroup, n_genes_genome, clade_groups, acc) %>%
#   distinct() %>%
#   group_by(clade_groups, Orthogroup) %>%
#   summarize(mean = median(n_genes_genome)) %>%
#   filter(clade_groups %in% c('Escovopsis', 'Trichoderma')) %>%
#   spread(clade_groups, mean) %>%
#   replace_na(list(Escovopsis = 0, Trichoderma = 0)) %>%
#   group_by(Escovopsis, Trichoderma) %>%
#   mutate(NumberGenes = length(Orthogroup)) %>%
#   ungroup() %>%
#   select(-Orthogroup) %>%
#   unique()
# 
# 
# line_data <- data.frame(y = seq(0, max(df$Escovopsis)), x = seq(0, max(df$Escovopsis)))
# 
# t_plot <- df %>%  ggplot(aes(x = Trichoderma, y = Escovopsis)) +
#   geom_line(data = line_data, aes(x = x, y=y)) +
#   geom_area(data = line_data, aes(x = x, y=y),
#             alpha = 0.5, position = 'identity') + 
#     geom_point(aes(size = NumberGenes)) +
#   theme_classic() +
#   labs(y = 'Median number of gene copies in Escovopsis',
#        x = 'Median number of gene copies in Trichoderma') +
#   scale_size(breaks = c(100, 1000, 3000), labels = c(100, 1000, 3000)) +
#   scale_x_continuous(breaks = seq(0, max(df$Trichoderma)), labels = seq(0, max(df$Trichoderma))) +
#   scale_y_continuous(breaks = seq(0, max(df$Escovopsis)), labels = seq(0, max(df$Escovopsis))) +
#   theme(
#     legend.position = c(0.05, .8), 
#     legend.justification = c(0, 0)
#   )
# t_plot
# 
# # with marginal histogram
# 
# 
# df <- orthog_meta_ls %>%
#   select(Orthogroup, n_genes_genome, clade_groups, acc) %>%
#   distinct() %>%
#   group_by(clade_groups, Orthogroup) %>%
#   summarize(mean = median(n_genes_genome)) %>%
#   filter(clade_groups %in% c('Escovopsis', 'Hypomyces_Cladobotryum')) %>%
#   spread(clade_groups, mean) %>%
#   replace_na(list(Escovopsis = 0, Hypomyces_Cladobotryum = 0)) %>%
#   group_by(Escovopsis, Hypomyces_Cladobotryum) %>%
#   mutate(NumberGenes = length(Orthogroup)) %>%
#   ungroup() %>%
#   select(-Orthogroup) %>%
#   unique()
# 
# line_data <- data.frame(y = seq(0, max(df$Escovopsis)), x = seq(0, max(df$Escovopsis)))
# 
# h_plot <- df %>%  ggplot(aes(x = Hypomyces_Cladobotryum, y = Escovopsis)) +
#   geom_line(data = line_data, aes(x = x, y=y)) +
#   geom_area(data = line_data, aes(x = x, y=y),
#             alpha = 0.5, position = 'identity') + 
#   geom_point(aes(size = NumberGenes)) +
#   theme_classic() +
#   ggtitle('') +
#   labs(y = 'Median number of gene copies in Escovopsis',
#        x = 'Median number of gene copies in Hypomyces_Cladobotryum') +
#   scale_size(breaks = c(100, 1000, 2000), labels = c(100, 1000, 2000)) +
#   scale_x_continuous(breaks = seq(0, max(df$Hypomyces_Cladobotryum)), labels = seq(0, max(df$Hypomyces_Cladobotryum))) +
#   scale_y_continuous(breaks = seq(0, max(df$Escovopsis)), labels = seq(0, max(df$Escovopsis))) +
#   theme(
#     legend.position = c(0.05, .8), 
#     legend.justification = c(0, 0)
#   )
# #h_plot
# 
# plot_grid(plotlist = list(t_plot, h_plot), labels = 'AUTO') %>%
#   ggsave2(filename = 'plots/copy_numer_scatter.pdf')



########################
##COMBINING OUTGROUPS##
#######################

# df <- orthog_meta_ls %>%
#   select(Orthogroup, n_genes_genome, clade_groups, acc, tool) %>%
#   distinct() %>%
#   mutate(clade_groups = ifelse(clade_groups %in% c('Trichoderma', 'Hypomyces_Cladobotryum'), 
#                                  yes = 'Outgroup',
#                                  no = clade_groups)) %>%
#   group_by(clade_groups, Orthogroup) %>%
#   summarize(mean = median(n_genes_genome)) %>%
#   spread(clade_groups, mean) %>%
#   replace_na(list(Escovopsis = 0, Outgroup = 0)) %>%
#   group_by(Escovopsis, Outgroup) %>%
#   mutate(NumberGenes = length(Orthogroup)) %>%
#   ungroup()
# 
# df  %>%
#   select(-Orthogroup) %>%
#   unique() %>%  
#   ggplot(aes(x = Outgroup, y = Escovopsis)) +
#   geom_line(data = line_data, aes(x = x, y=y), linetype = "dashed", colour = 'grey') +
# #  geom_area(data = line_data, aes(x = x, y=y),
# #            alpha = 0.5, position = 'identity') + 
#   geom_tile(aes(x = Outgroup, y = Escovopsis, fill = NumberGenes)) +
#   theme_classic() +
#   ggtitle('') +
#   labs(y = 'Median number of gene copies in Escovopsis',
#        x = 'Median number of gene copies in Outgroup') +
#   scale_size(breaks = c(100, 1000, 2000), labels = c(100, 1000, 2000)) +
#   scale_x_continuous(breaks = seq(0, max(df$Outgroup)), 
#                      labels = seq(0, max(df$Outgroup)), 
#                      expand = c(0, 0)) +
#   scale_y_continuous(breaks = seq(0, max(df$Escovopsis)), 
#                      labels = seq(0, max(df$Escovopsis)), 
#                      expand = c(0, 0)) 



df <- orthog_meta_ls %>%
  select(Orthogroup, n_genes_genome, clade_groups, acc, tool) %>%
  distinct() %>%
  mutate(clade_groups = ifelse(clade_groups %in% c('Trichoderma', 'Hypomyces_Cladobotryum'), 
                               yes = 'Outgroup',
                               no = clade_groups)) %>%
  group_by(clade_groups, tool, Orthogroup) %>%
  summarize(mean = median(n_genes_genome)) %>%
  spread(clade_groups, mean) %>%
  replace_na(list(Escovopsis = 0, Outgroup = 0)) %>%
  group_by(tool, Escovopsis, Outgroup) %>%
  mutate(NumberGenes = length(Orthogroup)) %>%
  ungroup() %>% 
  select(-Orthogroup) %>%
  distinct()

colors_heat <- colorRampPalette(c("#132B43", "#56B1F7"))

fig2a <- df %>%
  filter(tool != 'Single Copy Genes') %>%
  ggplot() +
  geom_point(aes(x = Outgroup, y = Escovopsis, size = NumberGenes)) +
  #geom_tile(aes(x = Outgroup, y = Escovopsis, fill = NumberGenes)) +
  geom_line(data = line_data, aes(x = x, y=y), linetype = "dashed", colour = 'grey') +
  facet_wrap(~tool, scales = 'free_x') +
  theme_classic() +
  labs(y = 'Median Copy Number in Escovopsis',
       x = 'Median Copy Number in Outgroups') +
  scale_x_continuous(breaks = seq(0, max(df$Outgroup), by = 2), 
                     labels = seq(0, max(df$Outgroup), by = 2)) +
  scale_y_continuous(breaks = seq(0, max(df$Escovopsis)), 
                     labels = seq(0, max(df$Escovopsis))) +
  # scale_color_gradientn(breaks = c(1, seq(100, 1600, by = 200)),
  #                      labels = c(1, seq(100, 1600, by = 200)),
  #                      colours = colors_heat(10),
  #                      values = c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,0.7, 0.8, 0.9, 1), 
  #                      guide="colourbar") +
  # guides(fill = guide_colourbar(direction = "horizontal", nbin = 50, title.vjust = 0.7, barwidth = 15)) +
  theme(legend.position = 'bottom')

fig2a

df %>%
  filter(tool == 'Single Copy Genes') %>%
  ggplot() +
  geom_tile(aes(x = Outgroup, y = Escovopsis, fill = NumberGenes)) +
  geom_line(data = line_data, aes(x = x, y=y), linetype = "dashed", colour = 'grey') +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, max(df$Outgroup)), 
                     labels = seq(0, max(df$Outgroup)), 
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, max(df$Escovopsis)), 
                     labels = seq(0, max(df$Escovopsis)), 
                     expand = c(0, 0)) +
  scale_fill_gradientn(breaks = c(1, seq(100, 1600, by = 200)),
                       labels = c(1, seq(100, 1600, by = 200)),
                       colours = colors_heat(10),
                       values = c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,0.7, 0.8, 0.9, 1), 
                       guide="colourbar") +
  guides(fill = guide_colourbar(barheight = 20, direction = "vertical",
                                title.position="top", title.hjust = 0.5,title.vjust = 0.5, nbin = 50))





##################################################################
#### TRYING DIVERSITY METRICS 
#################################################################

tools <- as.character(unique(orthog_meta_ls$tool))

div_mats <- lapply(tools, function(x){
  div_mat <- orthog_meta_ls %>% select(genome, tool, annotation) %>%
    filter(tool == x) %>%
    select(-tool) %>%
    group_by(genome, annotation) %>%
    mutate(n_annotation = length(annotation)) %>%
    ungroup() %>%
    distinct() %>%
    spread(annotation, n_annotation) %>%
    column_to_rownames(var = 'genome')
  div_mat[is.na(div_mat)] <- 0
  as.matrix(div_mat)
})

names(div_mats) <- tools



div_out <- lapply(names(div_mats), function(x){
  div_mat <- div_mats[[x]]
  H <- diversity(div_mat) ## shannon index
  evenness <- H/log(specnumber(div_mat)) ## evenness
  simp <- diversity(div_mat, "simpson")
  invsimp <- diversity(div_mat, "inv")
  S <- specnumber(div_mat) # observed number of species
  
  cbind(S, H, evenness, simp, invsimp) %>% 
    data.frame() %>%
    rownames_to_column(var = 'genome') %>%
    mutate(tool = x)
}) %>%
  bind_rows() %>%
  left_join(., metadata, by = 'genome')

div_out$Agriculture <- factor(div_out$Agriculture, levels = rev(c("Trichoderma",
                                                                                "Hypomyces_Cladobotryum",
                                                                                "Lower",
                                                                                "Coral",
                                                                                "Higher",
                                                                                "Leafcutter"
)))





          
div_out %>% ggplot(., aes(x = S, y = H, color = Agriculture)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors[na.omit(unique(metadata$Agriculture))]) +
  theme_classic() +
  labs(y = 'Entropy', x = 'Number of Genes') +
  facet_wrap(~tool, nrow = 3, ncol = 3, scales = 'free') +
  ggsave(filename = 'plots/shannon_diversity_gene_copy_scatter.pdf')

div_out %>%
  ggplot(., aes(x = Agriculture, y = H, color = Agriculture)) +
  geom_boxplot() +
  scale_color_manual(values = colors[na.omit(unique(metadata$Agriculture))]) +
  theme_classic() +
  facet_wrap(~tool, nrow = 3, ncol = 3, scales = 'free') +
  labs(y = 'Entropy', x = '') +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ggsave(filename = 'plots/shannon_diversity_gene_copy_box.pdf')


div_out %>%
  ggplot(., aes(x = Agriculture, y = S, color = Agriculture)) +
  geom_boxplot() +
  scale_color_manual(values = colors[na.omit(unique(metadata$Agriculture))]) +
  theme_classic() +
  facet_wrap(~tool, nrow = 3, ncol = 3, scales = 'free') +
  labs(y = 'Gene Richness', x = '') +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ggsave(filename = 'plots/richness_gene_copy_box.pdf')



div_out %>%
  ggplot(., aes(x = Agriculture, y = evenness, color = Agriculture)) +
  geom_boxplot() +
  scale_color_manual(values = colors[na.omit(unique(metadata$Agriculture))]) +
  theme_classic() +
  facet_wrap(~tool, nrow = 3, ncol = 3, scales = 'free') +
  labs(y = 'Gene Evenness', x = '') +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ggsave(filename = 'plots/evenness_gene_copy_box.pdf')



lapply(unique(div_out$tool), function(x){
  sub <- div_out %>% filter(tool == x)
# shan_p <- sub %>% ggplot(., aes(x = S, y = H, color = Agriculture)) +
#     geom_jitter(size = 2) +
#     scale_color_manual(values = colors[na.omit(unique(metadata$Agriculture))]) +
#     theme_classic() +
#   theme(legend.position = 'none') +
#     labs(y = 'Entropy', x = 'Number of Genes') 
  
  shan <-  sub %>%
    ggplot(., aes(x = Agriculture, y = H, color = Agriculture)) +
    geom_boxplot() +
    scale_color_manual(values = colors[na.omit(unique(metadata$Agriculture))]) +
    theme_classic() +
    theme(legend.position = 'none') +
    labs(y = 'Entropy', x = '') +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  rich <-  sub %>%
    ggplot(., aes(x = Agriculture, y = S, color = Agriculture)) +
    geom_boxplot() +
    scale_color_manual(values = colors[na.omit(unique(metadata$Agriculture))]) +
    theme_classic() +
    theme(legend.position = 'none') +
    labs(y = 'Gene Richness', x = '') +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  even <-  sub %>%
    ggplot(., aes(x = Agriculture, y = evenness, color = Agriculture)) +
    geom_boxplot() +
    scale_color_manual(values = colors[na.omit(unique(metadata$Agriculture))]) +
    theme_classic() +
    theme(legend.position = 'none') +
    labs(y = 'Gene Evenness', x = '') +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  all <- cowplot::plot_grid(plotlist = list(shan, rich, even), labels = 'AUTO')
  outfile = sub(" ", "-", paste0('plots/', x, '_gene_copy_box_diversity.pdf'))
  save_plot(plot = all, filename = outfile)
  
})


div_out$tool <- factor(div_out$tool, levels = c("Single Copy Genes",
                                                              "Resistance Genes",
                                                              "Virulence Genes",
                                                              "BGC",
                                                              "Lipase Genes",
                                                              "Peptidase Genes",
                                                              "CAZyme Genes"))

fig2b <- div_out %>% 
  filter(tool != 'Single Copy Genes') %>%
  ggplot(., aes(x = S, y = H, color = Agriculture)) +
  geom_jitter(size = 2) +
  scale_color_manual(values = colors[na.omit(unique(metadata$Agriculture))]) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  labs(y = 'Entropy', x = 'Number of Genes') +
  facet_wrap(~tool, scales = 'free')




cowplot::plot_grid(plotlist = list(fig2a, fig2b), labels = 'AUTO') %>%
  ggsave2( filename = 'plots/figure3.pdf')

### ONLY VIRULENCE CONVERGES< NOT CONSIDERING
# nmds_all <- lapply(tools, function(y){
#   x <- div_mats[[y]]
#   pca_out    <- prcomp(x)
#   eigs       <- pca_out$sdev^2
#   proportion <- (eigs/sum(eigs))*100
#   cumulative <- cumsum(eigs)/sum(eigs)
#   n_k        <- length(which(proportion >= 10))
#   
#   if(n_k < 2){
#     n_k = 2
#   }
#   print(paste(y, n_k))
#   
#   metaMDS(x, k = n_k) # The number of reduced dimensions ## components with >10% variance explained
# })
# names(nmds_all) <- tools
# 
# lapply(nmds_all, function(x){x$converged})
# ## make stressplots
# 
# lapply(tools, function(x){
#   stressplot(nmds_all[[x]])
# })
# 
# gene_pvals <- lapply(tools, function(y){
#   print(y)
#   nmds              <- nmds_all[[y]]
#   x                 <- div_mats[[y]]
#   vec.sp            <- envfit(nmds$points, x, perm=1000)
#   vec.sp.df         <- as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
#   vec.sp.df$pval    <- vec.sp$vectors$pval
#   vec.sp.df$annotation <- rownames(vec.sp.df) 
#   vec.sp.df$tool <- y
#   vec.sp.df
# })
# 
# 
# pvals_all <- gene_pvals %>% bind_rows()
# 
# write_csv(pvals_all, 'tables/nmds_significant_genes.csv')
# 
# pvals_all %>% filter(pval <= 0.001, tool != 'Virulence Genes', tool != 'Peptidase Genes', tool != 'Single Copy Genes') %>%
#   filter(tool == 'BGC') 
# 
# lapply(gene_pvals, function(x){
#   x %>% filter(pval <= 0.001) %>% nrow()
# })
# tools
# 
# y <- tools[[1]]
# nmds_ls <- lapply(tools, function(y){
#   nmds <- data.frame(scores(nmds_all[[y]])) %>%
#   rownames_to_column( var = 'genome') %>%
#   left_join(metadata, by = 'genome')
#   nmds %>%
#     ggplot() +
#     geom_point(aes(x = NMDS1, y = NMDS2, color = Agriculture), size = 2) +
#     scale_color_manual(values = colors[names(colors) %in% metadata$Agriculture]) +
#     theme_classic() +
#     theme(legend.position = 'none') +
#     labs(title = y) #+
#     # geom_segment(data=t,aes(x=0,xend=MDS1,y=0,yend=MDS2),
#     #              arrow = arrow(length = unit(0.3, "cm")),colour="grey",inherit_aes=FALSE) + 
#     # geom_text(data=t,aes(x=MDS1,y=MDS2,label=species),size=2)+
#     # coord_fixed()
# })
# nmds_ls[[1]]
# 
# pdf(file = 'plots/gene_copy_nmds.pdf')
# cowplot::plot_grid(plotlist = nmds_ls, align = 'hv')
# dev.off()

###############################################
### REPEAT ABOVE BUT WITHIN ANT AGRICULTURE ###
###############################################

ant_ag <- c('Coral', 'Higher', 'Lower', 'Leafcutter')

div_mats_ag <- lapply(tools, function(x){
  div_mat <- orthog_meta_ls %>% 
    filter(Agriculture %in% ant_ag) %>%
    select(genome, tool, annotation) %>%
    filter(tool == x) %>%
    select(-tool) %>%
    group_by(genome, annotation) %>%
    mutate(n_annotation = length(annotation)) %>%
    ungroup() %>%
    distinct() %>%
    spread(annotation, n_annotation) %>%
    column_to_rownames(var = 'genome')
  div_mat[is.na(div_mat)] <- 0
  as.matrix(div_mat)
})

names(div_mats_ag) <- tools



div_out <- lapply(names(div_mats_ag), function(x){
  div_mat <-div_mats_ag[[x]]
  H <- diversity(div_mat) ## shannon index
  evenness <- H/log(specnumber(div_mat)) ## evenness
  simp <- diversity(div_mat, "simpson")
  invsimp <- diversity(div_mat, "inv")
  S <- specnumber(div_mat) # observed number of species
  
  cbind(S, H, evenness, simp, invsimp) %>% 
    data.frame() %>%
    rownames_to_column(var = 'genome') %>%
    mutate(tool = x)
}) %>%
  bind_rows() %>%
  left_join(., metadata, by = 'genome')

div_out$Agriculture <- factor(div_out$Agriculture, levels = rev(c("Trichoderma",
                                                                  "Hypomyces_Cladobotryum",
                                                                  "Lower",
                                                                  "Coral",
                                                                  "Higher",
                                                                  "Leafcutter"
)))

nmds_all <- lapply(tools, function(y){
  x <-div_mats_ag[[y]]
  pca_out    <- prcomp(x)
  eigs       <- pca_out$sdev^2
  proportion <- (eigs/sum(eigs))*100
  cumulative <- cumsum(eigs)/sum(eigs)
  n_k        <- length(which(proportion >= 10))
  
  if(n_k < 2){
    n_k = 2
  }
  print(paste(y, n_k))
  
  mds <- metaMDS(x, k = n_k) # The number of reduced dimensions ## components with >10% variance explained
})

names(nmds_all) <- tools
lapply(nmds_all, function(x){x$converged})

## make stressplots

x <- tools[6]


lapply(tools, function(x){
  stressplot(nmds_all[[x]])
})


  

gene_pvals <- lapply(tools, function(y){
  print(y)
  nmds              <- nmds_all[[y]]
  x                 <-div_mats_ag[[y]]
  vec.sp            <- envfit(nmds$points, x, perm=1000)
  vec.sp.df         <- as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
  vec.sp.df$pval    <- vec.sp$vectors$pval
  vec.sp.df$annotation <- rownames(vec.sp.df) 
  vec.sp.df$tool <- y
  vec.sp.df
})
names(gene_pvals) <- tools

env <- metadata %>% select(genome, Agriculture1, coordinates, Host, cultivar) %>% na.omit() %>%
  filter(Agriculture %in% ant_ag) %>%
  column_to_rownames(var = 'genome')

y <- tools[5]

metadata_centroids <- lapply(tools, function(y){
  nmds              <- nmds_all[[y]]
  vec.sp            <- envfit(nmds$points, env, perm=1000)
  vec.sp.df <- cbind(rownames(vec.sp$factors$centroids), vec.sp$factors$centroids, vec.sp$factors$var.id) %>% data.frame(stringsAsFactors = FALSE)
  colnames(vec.sp.df) <- c('id', 'MDS1', 'MDS2', 'metadata_type')
  vec.sp.df$tool <- y
  vec.sp.df$pvals <- vec.sp$factors$pvals[vec.sp.df$metadata]
  vec.sp.df %>%
    mutate_at(c('MDS2', 'MDS1'), as.numeric)
}) %>% bind_rows()



all_counts <- do.call(cbind, div_mats) %>% t() %>%
  data.frame(., stringsAsFactors = FALSE) %>%
  rownames_to_column(var = 'annotation')

pvals_all <- gene_pvals %>% bind_rows() %>% left_join(all_counts, by = 'annotation')

write_csv(pvals_all, 'tables/nmds_significant_genes_ant_agriculture.csv')

pvals_all <- read_csv('tables/nmds_significant_genes_ant_agriculture.csv')

pvals_all %>% filter(pval <= 0.001)

lapply(gene_pvals, function(x){
  x %>% filter(pval <= 0.001) %>% nrow()
})





### MAKE HEATMAP OF COUNTS OF SIGNIFICANT GENES

pvals_all <- gene_pvals %>% bind_rows() %>% left_join(all_counts, by = 'annotation') %>%
  filter(pval <= 0.001, tool != 'Single Copy Genes') %>%
  select(-contains('MDS'), -pval)

col_annot <- pvals_all %>% select(annotation, tool) %>% column_to_rownames(var = 'annotation') %>%
  rename('Annotation' = tool)

col_annot$Annotation <- factor(col_annot$Annotation, levels = unique(col_annot$Annotation))

pvals_heat <- pvals_all %>% select(-tool) %>% column_to_rownames(var = 'annotation') %>% t()

pvals_heat_norm <- sweep(pvals_heat, 2, colSums(pvals_heat), `/`)



annot_colors        <- c('Virulence Genes' = '#ffc100', 
                         'Resistance Genes' = '#ff9a00',
                         'Lipase Genes' = '#1261A0', 
                         'CAZyme Genes' = '#3895D3', 
                         "Peptidase Genes" = '#58CCED')

annot_row <- metadata %>% select(genome, Agriculture) %>% na.omit() %>%
  column_to_rownames(var = 'genome')



#names(annot_colors) <- factor(unique(col_annot$Annotation), levels = unique(col_annot$Annotation))
annot_colors_ls     <- list(Annotation = annot_colors,
                            Agriculture = colors[unique(annot_row$Agriculture)])

result <- pvclust(pvals_heat_norm[order[-1], ], method.dist="binary",
                  method.hclust="average", nboot=100)

pheatmap(pvals_heat_norm[order[-1], ],
         color = colorRampPalette(brewer.pal(n = 5, name = "Greens"))(100),
         cluster_rows = FALSE,
         cluster_cols = result$hclust,
         show_colnames = F,
         annotation_col = col_annot,
         annotation_row = annot_row, 
         annotation_colors = annot_colors_ls,
         filename = 'plots/nmds_significant_genes_ant_agriculture_heat.pdf'
         ) 





ag_nmds_ls <- lapply(tools, function(y){
t <- gene_pvals[[y]] %>% filter(pval <= 0.001) %>%
    rename('genome' = annotation) %>%
    left_join(metadata, by = 'genome')
  nmds <- data.frame(scores(nmds_all[[y]])) %>%
    rownames_to_column( var = 'genome') %>%
    left_join(metadata, by = 'genome')
  p <- nmds %>%
    ggplot() +
    geom_point(aes(x = NMDS1, y = NMDS2, color = Agriculture), size = 2) +
    scale_color_manual(values = colors[names(colors) %in% metadata$Agriculture]) +
    theme_classic() +
    theme(legend.position = 'none') +
    labs(title = y)
  if(nrow(t) <= 20 && nrow(t) > 0){
    p+ geom_segment(data=t,aes(x=0, xend= MDS1, y=0, yend= MDS2),
                    arrow = arrow(length = unit(0.3, "cm")),colour="grey",inherit_aes=FALSE) +
      geom_text(data=t, aes(x=MDS1,y=MDS2,label=genome),size=2)
}else {p}
})

pdf(file = 'plots/gene_copy_ag_nmds.pdf')
cowplot::plot_grid(plotlist = ag_nmds_ls, align = 'hv')
dev.off()


ag_nmds_ls <- lapply(tools, function(y){
  t <- metadata_centroids %>% filter(tool == y, pvals <= 0.05, !grepl('Agriculture$', id), !grepl('Host', id)) %>%
    distinct()
  nmds <- data.frame(scores(nmds_all[[y]])) %>%
    rownames_to_column( var = 'genome') %>%
    left_join(metadata, by = 'genome')
  p <- nmds %>%
    ggplot() +
    geom_point(aes(x = NMDS1, y = NMDS2, color = Agriculture), size = 2) +
    scale_color_manual(values = colors[names(colors) %in% metadata$Agriculture]) +
    theme_classic() +
    theme(legend.position = 'none') +
    labs(title = y)
  if(nrow(t) <= 20 && nrow(t) > 0){
    p+ geom_segment(data=t,aes(x=0, xend= MDS1, y=0, yend= MDS2),
                    arrow = arrow(length = unit(0.3, "cm")),colour="grey",inherit_aes=FALSE) +
      geom_text(data=t, aes(x=MDS1,y=MDS2,label=id),size=2)
  }else {p}
})

pdf(file = 'plots/gene_copy_ag_nmds_metadata_centroid.pdf')
cowplot::plot_grid(plotlist = ag_nmds_ls, align = 'hv')
dev.off()
###############################################################
###  Trying plotting genes, largely a fail except for BGC  ####
###############################################################
# gene_pvals <- lapply(tools, function(y){
#   print(y)
#   nmds              <- nmds_all[[y]]
#   x                 <- t(div_mats_ag[[y]])
#   vec.sp            <- envfit(scores(nmds, display = 'species'), x, perm=1000)
#   vec.sp.df         <- as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
#   vec.sp.df$pval    <- vec.sp$vectors$pval
#   vec.sp.df$annotation <- rownames(vec.sp.df)
#   vec.sp.df$tool <- y
#   vec.sp.df
# })
# names(gene_pvals) <- tools


y = tools[7]
lapply(tools, function(y){
  t <- metadata_centroids %>% filter(tool == y, pvals <= 0.05, grepl('Agriculture1', id)) %>%
    mutate(id = sub('Agriculture1', '', id))
  nmds <- data.frame(scores(nmds_all[[y]], display = 'species')) %>%
    rownames_to_column( var = 'annotation') %>%
    left_join(annot, by = 'annotation') %>%
    left_join(metadata, by = 'genome') %>%
    distinct() %>%
    filter(Agriculture %in% ant_ag) %>%
    select(-gene) %>%
    distinct()
  p <- nmds %>%
    ggplot(aes(x = NMDS1, y = NMDS2)) +
    geom_jitter(height = 0.1, width = 0.12, size = 2, alpha = 0.2, aes(color = Agriculture)) +
    scale_color_manual(values = colors[names(colors) %in% metadata$Agriculture]) +
    theme_classic() +
    labs(title = paste(y, 'Genes'))
    if(nrow(t) > 0){
    p+ geom_segment(data=t,aes(x=0, xend= MDS1, y=0, yend= MDS2),
                     arrow = arrow(length = unit(0.3, "cm")),colour="grey",inherit_aes=FALSE) +
      geom_text(data=t, aes(x=MDS1,y=MDS2,label=id),size=2)+
      coord_fixed()
    }else {p}
})





