library(tidyverse)
library(pheatmap)
library(vegan)
library(ggord)
library(dendextend)

munge_gca <- function(x){
  paste(strsplit(x, split = "_")[[1]][c(1,2)], collapse = "_")
}


colors <- c("Coral" = "#CE3DD0",
                "Higher" = "#2D71F6",
                "Lower" = "#b8860b",
                "Leafcutter" = "#377D22",
                "Outgroup1" = 'black',
                "Outgroup2" = 'gray')


metadata <- read_csv('tables/metadata.csv')
meta_levels <- metadata$genus_species

## dive in deeper to presence absence output
ant_ag <- c("Lower",
            "Coral",
            "Higher", 
            "Leafcutter")



files <- list.files("bigscape/parsed_networks", 
                    pattern = 'network', 
                    recursive = TRUE, 
                    full.names = TRUE)

all_data <- lapply(files, function(x){
  table <- read_tsv(x, col_names = c('component', 'bigscape')) %>%
    mutate(BGC_type = sub("_.*$", "", basename(x)))
}) %>% bind_rows() %>%
  group_by(bigscape) %>%
  mutate(acc = case_when(
    grepl('GCA', bigscape) ~ munge_gca(bigscape),
    grepl('SPDT', bigscape) ~ "SPDT00000000",
    TRUE ~ sub("_.*$", "", bigscape)
  )) %>%
  ungroup() %>%
  left_join(., metadata, by = 'acc') %>%
  mutate(Presence = 1) %>%
  group_by(component, BGC_type) %>%
  mutate(num_agricultures = length(unique(Agriculture[which(Agriculture %in% ant_ag)]))) %>% ## removes outgroups from count
  filter(!is.na(acc)) %>%
  mutate(Agriculture = ifelse(is.na(Agriculture),
                              yes = 'Outgroup',
                              no = Agriculture)) %>%
  ungroup() %>%
  group_by(BGC_type, component) %>%
  mutate(bgc_plot_label = ifelse(any(grepl('BGC', acc)),
                                 yes = paste(grep('BGC', acc, value = TRUE), collapse = ","),
                                 no = paste0(BGC_type, "_", component))) %>%
  ungroup()


#View(select(all_data, component, Agriculture, bigscape, num_agricultures, BGC_type))


#unique(all_data$acc)

#all_data %>% filter(component == 'component_111', BGC_type == 'mix') %>% View()

## do all bgc first with 'mix' seperate


ann_colors = list(
  'BGC_type' = c(NRPS = '#1B9E77', PKSI = '#D95F02', Terpene = '#7570B3', Others = "#E7298A", PKS.NRP = "#66A61E"),
  'Agriculture' = c(Lower = '#FFFEAB',
                      Coral = "#CE3DD0",
                      Higher = "#2D71F6",
                      Leafcutter = "#377D22", Outgroup2 = 'gray', Outgroup1 = 'black'))


all_bgc <- all_data %>%
  filter(!BGC_type %in% c('PKSother', 'mix'),
         num_agricultures > 0,
         !grepl('BGC', acc)) %>% ## only look at BGC that are in at least 1 ant agriculture
  select(bgc_plot_label, genus_species, Presence, BGC_type) %>%
  mutate(bgc_plot_label = sub('-', '.', bgc_plot_label)) %>%
  mutate(BGC_type = sub('-', '.', BGC_type))



all_data %>% filter(num_agricultures >= 3) %>%
  distinct() %>%
  rowwise() %>%
  mutate(file = grep(paste0(component, ".gff3"), list.files('bigscape/gff3', pattern = BGC_type, full.names = TRUE), value = TRUE)) %>%
  select(file) %>%
  distinct() %>%
  write_tsv(., path = 'tables/bgc_greater_than_2_ags.txt', col_names = FALSE)

all_bgc_s <- all_bgc %>% select(-BGC_type) %>%
  distinct() %>%
  spread(bgc_plot_label, Presence) %>%
  data.frame()

rownames(all_bgc_s) <- all_bgc_s$genus_species

all_bgc_s$genus_species <- NULL
all_bgc_s <- as.matrix(all_bgc_s)
all_bgc_s[is.na(all_bgc_s)] <- 0

## rows are genomes
rows_order <- meta_levels[which(meta_levels %in% rownames(all_bgc_s))]

m <- metadata %>% select(Agriculture, genus_species) %>%  mutate(Agriculture = replace_na(Agriculture, 'Outgroup')) %>%
  distinct() %>% data.frame()
rownames(m) <- m$genus_species
m$genus_species <- NULL




annotation_col <- all_bgc %>%
  select(-genus_species, -Presence) %>%
  filter(bgc_plot_label %in% colnames(all_bgc_s)) %>%
  distinct() %>%
  data.frame()

rownames(annotation_col) <- annotation_col$bgc_plot_label
annotation_col$bgc_plot_label <- NULL


annotation_row = data.frame(
     "Agriculture" = m[rownames(all_bgc_s), ],
     row.names = rownames(all_bgc_s),
    stringsAsFactors = FALSE
   )


pheatmap::pheatmap(all_bgc_s[rows_order,],
                   legend = FALSE,
                   color = c('grey87', 'black'),
                   cluster_rows = clusters,
                   cluster_cols = TRUE,
                   border_color = "grey70",
                   cellwidth = 10,
                   cellheight = 10,
                   annotation_row = annotation_row,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   filename = paste0('plots/all_BGC_components_heatmap.pdf'))



## Ordinate all BGCs with no threshold cutoff

all_bgc <- all_data %>%
  filter(!BGC_type %in% c('mix'),
         !grepl('BGC', acc)) %>% ## only look at BGC that are in at least 1 ant agriculture
  select(bgc_plot_label, genus_species, Presence, BGC_type) %>%
  mutate(bgc_plot_label = sub('-', '.', bgc_plot_label)) %>%
  mutate(BGC_type = sub('-', '.', BGC_type))

#all_bgc %>% filter(component == 'component_111', BGC_type == 'mix') %>% View()

all_bgc_s <- all_bgc %>% select(-BGC_type) %>%
  distinct() %>%
  spread(bgc_plot_label, Presence) %>%
  data.frame()


rownames(all_bgc_s) <- all_bgc_s$genus_species
all_bgc_s$genus_species <- NULL
all_bgc_s <- as.matrix(all_bgc_s)
all_bgc_s[is.na(all_bgc_s)] <- 0




library(pvclust)
set.seed(1234)
result <- pvclust(t(all_bgc_s), method.dist="cor",
                  method.hclust="average", nboot=100)

pdf(file = 'plots/bootstrapped_BGC_dendrogram.pdf')
plot(result)
pvrect(result, pv = 'bp')
dev.off()

dend_colors        <- metadata$Agriculture
names(dend_colors) <- metadata$genus_species


hcd <- as.dendrogram(result)  %>%
  set("leaves_pch", 19) %>%
  set("leaves_col", colors[dend_colors[result$hclust$labels[result$hclust$order]]])



pdf(file = 'plots/bootstrapped_BGC_dendrogram_colors.pdf')
par(mar = c(2, 2, 2, 10))
plot(hcd, horiz = TRUE, edgePar = list(lwd = 2))
dev.off()



pca_data <- prcomp(all_bgc_s)
eigs <- pca_data$sdev^2
proportion = (eigs/sum(eigs))*100
head(proportion)
cumulative = cumsum(eigs)/sum(eigs)
screeplot(pca_data)


n_k <- length(which(proportion >= 10))

## nmds
## bray curtis distance is more resilient to nulls
example_NMDS=metaMDS(all_bgc_s, k = 2) # The number of reduced dimensions ## components with >10% variance explained

stressplot(example_NMDS)


treat = replace_na(annotation_row[,'Agriculture'], replace = 'Outgroup') ## corresponds to rows/communities (genomes)


## test for significance
bgc_dist <- vegdist(all_bgc_s)
ano_test <- anosim(bgc_dist, grouping = treat)
summary(ano_test)
plot(ano_test)

# ggord(example_NMDS, grp_in = treat, arrow = NULL, obslab = FALSE,
#       txt = FALSE,
#       poly=FALSE, size=2) + theme_classic() +
#   labs(title = 'BGCs stratify by ant agriculture',
#        subtitle = paste('Pval:', ano_test$signif, ", R:", round(ano_test$statistic, digits = 2))) +
#   ggsave('plots/BGC_all_ord.pdf')



data2 <- example_NMDS$points %>% 
  data.frame() %>% 
  rownames_to_column(var = 'genus_species') %>%
  left_join(., metadata, by = 'genus_species')

ggplot(data2, aes(x = MDS1, y = MDS2)) + 
  geom_point(shape = 21, size = 3, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  theme_bw() +
  geom_text(data = filter(data2, genus_species %in% c('ICBG712', 'ICBG721')), 
            aes(label = genus_species, color = NULL), nudge_x = -0.04, hjust = 0) + 
  scale_fill_manual(values = colors) +
  ggsave('plots/BGC_all_ord.pdf', width = 7, height = 5)





## Make the same plots for the 'mix' analysis



all_bgc <- all_data %>%
  filter(BGC_type %in% c('mix'),
         num_agricultures > 0,
         !grepl('BGC', acc)) %>% ## only look at BGC that are in at least 1 ant agriculture
  select(bgc_plot_label, genus_species, Presence, BGC_type) %>%
  mutate(bgc_plot_label = sub('-', '.', bgc_plot_label))

#all_bgc %>% filter(component == 'component_111', BGC_type == 'mix') %>% View()

all_bgc_s <- all_bgc %>% select(-BGC_type) %>%
  distinct() %>%
  spread(bgc_plot_label, Presence) %>%
  data.frame()

rownames(all_bgc_s) <- all_bgc_s$genus_species
all_bgc_s$genus_species <- NULL
all_bgc_s <- as.matrix(all_bgc_s)
all_bgc_s[is.na(all_bgc_s)] <- 0

rows_order <- meta_levels[which(meta_levels %in% rownames(all_bgc_s))]

m <- metadata %>% select(Agriculture, genus_species) %>% mutate(Agriculture = replace_na(Agriculture, 'Outgroup')) %>%
  distinct() %>% data.frame()
rownames(m) <- m$genus_species
m$genus_species <- NULL




annotation_row = data.frame(
  "Agriculture" = m[rownames(all_bgc_s), ],
  row.names = rownames(all_bgc_s)
)

pheatmap::pheatmap(all_bgc_s[rows_order,],
                   legend = FALSE,
                   color = c('grey87', 'black'),
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   border_color = "grey70",
                   cellwidth = 10,
                   cellheight = 10,
                   annotation_row = annotation_row,
                   annotation_colors = ann_colors,
                   filename = paste0('plots/all_BGC_components_mix_heatmap.pdf'))
## Ordinate all BGCs

all_bgc <- all_data %>%
  filter(BGC_type %in% c('mix'),
         !grepl('BGC', acc)) %>% ## only look at BGC that are in at least 1 ant agriculture
  select(bgc_plot_label, genus_species, Presence, BGC_type) %>%
  mutate(bgc_plot_label = sub('-', '.', bgc_plot_label)) %>%
  mutate(BGC_type = sub('-', '.', BGC_type))

#all_bgc %>% filter(component == 'component_111', BGC_type == 'mix') %>% View()

all_bgc_s <- all_bgc %>% select(-BGC_type) %>%
  distinct() %>%
  spread(bgc_plot_label, Presence) %>%
  data.frame()

rownames(all_bgc_s) <- all_bgc_s$genus_species
all_bgc_s$genus_species <- NULL
all_bgc_s <- as.matrix(all_bgc_s)
all_bgc_s[is.na(all_bgc_s)] <- 0


pca_data <- prcomp(all_bgc_s)


eigs <- pca_data$sdev^2
proportion <- (eigs/sum(eigs))*100
head(proportion)
cumulative <- cumsum(eigs)/sum(eigs)
screeplot(pca_data)


n_k <- length(which(proportion >= 10))

## nmds
## bray curtis distance is more resilient to nulls
example_NMDS=metaMDS(all_bgc_s, k = n_k) # The number of reduced dimensions ## components with >10% variance explained

stressplot(example_NMDS)


treat = replace_na(annotation_row[,'Agriculture'], replace = 'Outgroup') ## corresponds to rows/communities (genomes)


## test for significance
bgc_dist <- vegdist(all_bgc_s)
ano_test <- anosim(bgc_dist, grouping = treat)
summary(ano_test)
plot(ano_test)

# ggord(example_NMDS, grp_in = treat, arrow = NULL, obslab = FALSE,
#       txt = FALSE,
#       poly=FALSE, size=2) + theme_classic() +
#   labs(title = 'BGCs stratify by ant agriculture',
#        subtitle = paste('Pval:', ano_test$signif, ", R:", round(ano_test$statistic, digits = 2))) +
#   ggsave('plots/BGC_mix_ord.pdf')

data2 <- example_NMDS$points %>% 
  data.frame() %>% 
  rownames_to_column(var = 'genus_species') %>%
  left_join(., metadata, by = 'genus_species')

ggplot(data2, aes(x = MDS1, y = MDS2)) + 
  geom_point(shape = 21, size = 3, aes(fill = Agriculture), stroke = 0.5, colour = 'black') +
  theme_bw() +
  geom_text(data = filter(data2, genus_species %in% c('ICBG712', 'ICBG721')), 
            aes(label = genus_species, color = NULL), nudge_x = -0.04, hjust = 0) + 
  scale_fill_manual(values = colors) +
  ggsave('plots/BGC_mix_ord.pdf', width = 7, height = 7)





