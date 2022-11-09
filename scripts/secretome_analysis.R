library(tidyverse)
library(vegan)
library(RColorBrewer)
library(pheatmap)
library(pvclust)
#library("googlesheets4")


## just look at cazy num genes v num domains plot
# ag_colours <- c("Generalist" = "#808080",
#                 "Coral" = "#CE3DD0",
#                 "Higher" = "#2D71F6",
#                 "Lower" = "#FFFEAB",
#                 "Leafcutter" = "#377D22",
#                 "Outgroup1" = 'black',
#                 "Outgroup2" = 'gray')


shorten_genus_species <- function(x){
  spl = strsplit(x, split = " ")
  paste0(strsplit(spl[[1]][1], split = "")[[1]][1], '. ', spl[[1]][2])
}

## lipases
led_f <- list.files('annotation/LED', full.names = TRUE)



#paste(strsplit(basename(x), split = "_")[[1]][c(1,2)], collapse = "_")



led <- lapply(led_f, function(x){
  read_tsv(x, skip = 4, col_names = c('gene', 'lipase')) %>%
    mutate(genome_f = basename(x)) %>%
    mutate(genome = paste(strsplit(genome_f, split = "_")[[1]][c(1,2)], collapse = "_")) %>%
    mutate(genome = sub("\\..*$", "", genome))

}) %>% bind_rows() %>%
  select(gene, lipase, genome) %>%
  rename('value' = lipase) %>%
  mutate(type = 'led') %>%
  mutate(value = sub('\\..*$', '', value)) %>% distinct()


# led %>% group_by(genome_f) %>%
#   summarize(length(unique(lipase)),
#             length(unique(gene))) %>%
#   View()

## peptidase

merops_f <- list.files('annotation/merops', full.names = TRUE)

merops <- lapply(merops_f, function(x){
  read_tsv(x, col_names = c('gene', 'peptidase', 'eval', 'p_id', 'len')) %>%
    mutate(genome_f = basename(x))    %>%
    mutate(genome = paste(strsplit(genome_f, split = "_")[[1]][c(1,2)], collapse = "_")) %>%
    mutate(genome = sub("\\..*$", "", genome)) %>%
    mutate(peptidase = sub("-.*$", "", peptidase))
}) %>% bind_rows() %>%
  select(gene, peptidase, genome) %>%
  rename('value' = peptidase) %>%
  mutate(type = 'merops')


viru_f <- list.files('annotation/virulence/parsed', full.names = TRUE, pattern = 'fa')

viru <- lapply(viru_f, function(x){
  read_tsv(x, col_names = c('gene', 'virulence', 'eval', 'p_id', 'len')) %>%
    mutate(genome_f = basename(x))    %>%
    mutate(genome = paste(strsplit(genome_f, split = "_")[[1]][c(1,2)], collapse = "_")) %>%
    mutate(genome = sub("\\..*$", "", genome))
}) %>% bind_rows() %>%
  select(gene, virulence, genome)  %>%
  rename('value' = virulence) %>%
  mutate(type = 'virulence')


## cazymes



cazy_files <- list.files('annotation/cazy', full.names = TRUE)

#c <- cazy_files[1]

cazy <- lapply(cazy_files, function(c){
  read_tsv(c, col_names = c('gene', 'HMMER', 'Hotpep', 'DIAMOND', 'Signalp', '#ofTools'), skip = 1, col_types = 'cccccn') %>%
    mutate(genome = sub(".txt", "", basename(c))) %>%
    mutate(cazy = 1) %>%
    mutate(HMMER = sub("\\(.*$", '', HMMER))
}) %>% bind_rows() %>%
  distinct() %>%
  group_by(genome, HMMER) %>%
  mutate(n_genes = as.numeric(length(unique(gene)))) %>%
  ungroup()  %>%
  rename('value' = HMMER) %>%
  mutate(type = 'cazy') %>%
  mutate(genome = sub("\\..*$", "", genome)) %>%
  filter(value != "-")



# interpro_f <- list.files('annotation/interpro/', full.names = TRUE)
# 
# interpro <- lapply(interpro_f, function(x){
#   read_tsv(x) %>%
#     mutate(genome_f = basename(x))    %>%
#     mutate(genome = paste(strsplit(genome_f, split = "_")[[1]][c(1,2)], collapse = "_")) %>%
#     mutate(genome = sub("\\..*$", "", genome))
# }) %>% bind_rows()
# 
# ## filter IDs that are likely secreted proteins
# interpro_filtered <- interpro %>%
#   select(genome, Gene, SignalP_EUK, TMHMM) %>%
#   filter(!is.na(TMHMM) | !is.na(SignalP_EUK)) %>%
#   gather(tool, value, -genome, -Gene) %>%
#   filter(!is.na(value)) %>%
#   group_by(genome, Gene, value) %>%
#   mutate(count = length(value)) %>%
#   select(-tool) %>%
#   spread(value, count) %>%
#   rename('gene' = Gene) %>%


targetp_f <- list.files('annotation/targetp', full.names = TRUE)

# Prediction of localization, based on the scores above; the possible values are:
#   C	Chloroplast, i.e. the sequence contains cTP, a chloroplast transit peptide;
# M	Mitochondrion, i.e. the sequence contains mTP, a mitochondrial targeting peptide;
# S	Secretory pathway, i.e. the sequence contains SP, a signal peptide;


targetp <- lapply(targetp_f, function(x){
  read_tsv(x, comment = '#', col_names = c('gene', 'prediction', 'noTP')) %>%
    mutate(genome_f = basename(x))    %>%
    mutate(genome = paste(strsplit(genome_f, split = "_")[[1]][c(1,2)], collapse = "_")) %>%
    mutate(genome = sub("\\..*$", "", genome))
}) %>% bind_rows() %>%
  filter(prediction == 'SP') %>% ## pull out secreted proteins
  select(gene, genome)


card_f <- list.files('annotation/card', full.names = TRUE)


card <- lapply(card_f, function(x){
  read_tsv(x)  %>% select(ORF_ID, ARO) %>%
    rename('gene' = ORF_ID) %>%
    rename('value' = ARO) %>%
    mutate(value = as.character(value)) %>%
    mutate(type = 'card') %>%
    mutate(gene = sub(" .*$", '', gene)) %>%
    mutate(genome = basename(x)) %>%
    mutate(genome = case_when(
      grepl('GCA', genome) ~ sub("\\..*$", "", genome),
      grepl("LGSR", genome) ~ sub('.proteins.fasta.txt', '', genome),
      grepl('SPDT', genome) ~ sub('\\.1_genomic.all.maker.proteins.fasta.txt', '', genome),
      TRUE ~ sub('.all.maker.proteins.fasta.txt', '', genome)
    ))
}) %>% bind_rows()  %>%
  mutate(genome = sub("\\..*$", "", genome))


#wolfp_f <- list.files('wolfpsort', full.names = TRUE)

# cazy_spread <- cazy %>%
#   select(gene, genome, cazyme)

#cazy_spread[is.na(cazy_spread)] <- 0
all_annotations <- read_tsv('annotation/all_annotations.txt', col_names = c('genome', 'gene', 'tool', 'value'))

## number of annotations
all_annotations %>% select(genome, gene) %>% distinct() %>% .$genome %>% table()

interpro_filtered <- all_annotations %>% 
  filter(tool %in% c('TMHMM', 'SignalP_EUK')) %>% 
  tidyr::spread(key = tool, value = value) %>% 
  filter(!is.na(SignalP_EUK) | TMHMM == 1) %>%
  .$gene


## I'm looking for extr

## cytoskeleton, cytoplasm, nucleus, mitochondria, vesicles of secretory system, 
## endoplasmic reticulum (ER), Golgi, vacuole, plasma membrane, peroxisome, 
## extracellular space including cell wall.

wolfpsort <- all_annotations %>% filter(tool == 'wolfpsort', grepl('extr', value)) %>%
  select(gene, genome) %>%
  mutate(genome = sub("\\..*$", "", genome))



all_secreted <- bind_rows(targetp, wolfpsort) %>% distinct()




all_enzymes <- cazy %>%
  select(gene, genome, type, value) %>%
  filter(value != 'N') %>%
  bind_rows(., merops, led, viru, card) %>%
#  gather(type, value, -gene, -genome) %>% 
  filter(!is.na(value)) %>%
  mutate(count = 1)


secreted <- all_secreted %>%
  filter(gene %in% interpro_filtered) %>%
  left_join(., all_enzymes, by = c('genome', 'gene')) %>%
  filter(!is.na(type)) %>%
  mutate(genome = ifelse(genome == 'SPDT00000000', yes = 'GCA_008477525', no = genome))

n_genomes_total <- n_distinct(secreted$genome)

### %% percentage shared vs. not shared for each of the big groups


metadata <- read_csv('tables/metadata.csv') %>% 
  mutate(genome = sub("\\..*$", "", acc))




sum1 <- secreted %>% 
  group_by(value) %>% 
  mutate(genome_perc = length(unique(genome))/n_genomes_total) %>%
  left_join(metadata, by = 'genome')





sum1 %>% select(genome_perc, type) %>%
  distinct() %>%
  ggplot(aes(x = type, y = genome_perc)) + geom_boxplot(notch = TRUE) + theme_classic() +
  labs(y = 'Percentage of Genomes with a Given Annotation', x = 'Categories of Annotations')


sum1 %>% ggplot(aes(x = type, y = genome_perc)) + geom_boxplot(notch = TRUE) + theme_classic() +
  labs(y = 'Percentage of Genomes with a Given Annotation', x = 'Categories of Annotations')

sum1 %>%
  group_by(type) %>%
  summarize(length(unique(value)))

# can you also distinguish what % is in 95% of Esco but not others and 95% of other but not Esco for each of these?

n_genomes <- n_distinct(secreted$genome)

gene_representation_df <- secreted %>% 
  left_join(metadata, by = 'genome') %>%
  group_by(type, value) %>%
  mutate(representation_all = length(unique(genome))/n_genomes) %>% 
  group_by(Agriculture) %>%
  mutate(n_genomes_ag= length(unique(genome))) %>%
  group_by(Agriculture, type, value) %>%
  mutate(representation_ag = length(unique(genome))/n_genomes_ag) %>% 
  ungroup() %>%
  distinct() %>%
  select(Agriculture, value, contains('representation')) %>% distinct() %>%
  rename('id' = value) 


gene_representation_df %>% 
  mutate(genus = case_when(
    Agriculture == 'Trichoderma' ~ Agriculture,
    Agriculture == 'Hypomyces_Cladobotryum' ~ Agriculture,
    TRUE ~ 'Escovopsis'
  )) %>%
  group_by(genus) %>%
  summarize(length(unique(id)))


gene_representation_df <- gene_representation_df %>% 
  spread(Agriculture, representation_ag)

is.na(gene_representation_df) <- 0







# df %>% ggplot(aes(x = type, y = representation)) + geom_boxplot(notch = TRUE) + theme_classic() +
#   labs(y = 'Percentage of Genomes with a Given Annotation', x = 'Categories of Annotations') +
#   facet_wrap(~genus_group)
  

secreted %>% left_join(metadata, by = 'genome') %>%
  filter(value == '3003392') %>% View()



  

secreted_metadata <-  secreted %>% left_join(metadata, by = 'genome')  %>%
  group_by(genus) %>% 
  mutate(n_genomes = length(unique(genome))) %>%
  group_by(genus, value) %>% 
  mutate(genome_frac = length(unique(genome))/n_genomes) %>%
  ungroup()


## present in Trichoderma

#secreted_metadata %>%
#  filter(genus == 'Trichoderma') %>% .$value %>% n_distinct()

secreted_metadata %>%
  filter(genus == 'Trichoderma', genome_frac >= .90) %>% .$value %>% n_distinct()
 
## present in Escovopsis
#secreted_metadata %>%
#  filter(genus == 'Escovopsis') %>% .$value %>% n_distinct()

secreted_metadata %>%
  filter(genus == 'Escovopsis', genome_frac >= .90) %>% .$value %>% n_distinct()

## present in Hypomyces
#secreted_metadata %>%
#  filter(genus == 'Hypomyces') %>% .$value %>% n_distinct()

secreted_metadata %>%
  filter(genus == 'Hypomyces', genome_frac >= .90) %>% .$value %>% n_distinct()

## present in Cladobotryum

#secreted_metadata %>%
#  filter(genus == 'Cladobotryum') %>% .$value %>% n_distinct()

secreted_metadata %>%
  filter(genus == 'Cladobotryum', genome_frac >= .90) %>% .$value %>% n_distinct()




###########################################################################################
######################### CLUSTERING ANALYSIS #############################################
###########################################################################################

secreted_wide <- secreted %>%
  group_by(genome, type, value, gene) %>%
  mutate(count = sum(count)) %>%
  ungroup() %>%
  select(-gene) %>%
  distinct() %>%
  select(value, genome, count) %>% 
  spread(value, count) %>%
  column_to_rownames('genome')

secreted_wide[is.na(secreted_wide)] = 0

secreted_wide <- secreted_wide[,which(colSums(secreted_wide) > 1)]

#secreted_wide[is.na(secreted_wide)] = 0




bgc <- all_annotations %>% filter(tool == 'antismash') %>%
  filter(gene %in% secreted$gene)



col_annot <- secreted %>% 
  # mutate(in_BGC = ifelse(gene %in% bgc$gene, yes = TRUE, no = FALSE)) %>%
  # group_by(value) %>%
  # mutate(num_in_BGC = length(which(in_BGC)),
  #        num_total = length(in_BGC)) %>%
  # mutate(ratio = num_in_BGC*100/num_total) %>%
  # ungroup()
  select(value, type) %>%
  distinct() %>%
  filter(value %in% colnames(secreted_wide)) %>%
  column_to_rownames(var = 'value') %>%
  rename('Annotation' = type)






col_annot$Annotation <- factor(col_annot$Annotation, levels = unique(col_annot$Annotation))



annot_colors = brewer.pal(length(unique(col_annot$Annotation)), 'Dark2')
names(annot_colors)        <- unique(col_annot$Annotation)




annot_row <- metadata %>% select(genome, Agriculture) %>%
  na.omit() %>%
  column_to_rownames(var = 'genome')





#names(annot_colors) <- factor(unique(col_annot$Annotation), levels = unique(col_annot$Annotation))
source('scripts/color_palettes.R')




#result_col <- pvclust(secreted_wide, method.dist="binary",
#                 method.hclust="ward.D", nboot=1000)


#result_row <- pvclust(t(secreted_wide), method.dist="binary",
#                      method.hclust="ward.D", nboot=1000)



load('rdata/secretome.Rdata')

pdf('plots/secretome_col_pvclust.pdf', width = 60)
plot(result_col, labels = FALSE)
lines(result_col)
dev.off()




pdf('plots/secretome_row_pvclust.pdf', width = 10)
plot(result_row, labels = FALSE)
lines(result_row)
dev.off()





significant_cols <- pvpick(result_col)$clusters %>% unlist()
significant_rows <- pvpick(result_row)$clusters %>% unlist()

annot_colors_ls     <- list(Annotation = annot_colors,
                            Agriculture = colors[unique(annot_row$Agriculture)],
                            Significant = c('TRUE' = 'black', 'FALSE' = 'white'))







#annot_row$Significant <- ifelse(rownames(annot_row) %in% significant_rows, 'TRUE', 'FALSE')
#col_annot$Significant <- ifelse(rownames(col_annot) %in% significant_cols, 'TRUE', 'FALSE')

library(dendextend)


order_row <- c("GCA_003025115", "GCA_003025155", "GCA_001050175", "GCA_000167675", "GCA_000513815", "GCA_001481775", 
  "GCA_002894205", "GCA_000171015", "GCA_003025105", "GCA_002894145", "GCA_002838845", "GCA_003025095", 
  "GCA_000988865", "GCA_002022785", "GCA_011066345", "GCA_000170995", "GCA_003012105",   "GCA_008477525", "GCA_004303015", "GCA_011799845", "ICBG712", "ICBG721", 
 "ICBG2047", "ICBG2049", "ICBG2046", 
  "ICBG2048", "NIGD00000000", "ICBG710", "ICBG730", "ICBG1054", "ICBG726", "ICBG1065", "ICBG1075",
 "ICBG751", 
 "ICBG1096", "ICBG733", "ICBG742", "NIGB00000000", "LGSR00000000", "NQYS00000000", "ICBG731", 
 "ICBG736", "NQYQ00000000", "NIGC00000000", "NQYR00000000"
)



rotated_row <- dendextend::rotate(result_row$hclust, order = order_row)




pheatmap(secreted_wide,
         color = colorRampPalette(brewer.pal(n = 5, name = "Greens"))(100),
         cluster_rows = rotated_row,
         cluster_cols = result_col$hclust,
         show_colnames = F,
         annotation_col = col_annot,
         annotation_row = annot_row, 
         annotation_colors = annot_colors_ls,
         cellheight = 9,
         cellwidth = 1,
         filename = 'plots/secretome_heatmap.pdf'
) 


pheatmap(t(secreted_wide),
         color = colorRampPalette(brewer.pal(n = 5, name = "Greens"))(100),
         cluster_cols = rotated_row,
         cluster_row = result_col$hclust,
         #cluster_cols = result$hclust,
         annotation_row = col_annot,
         annotation_col = annot_row, 
         annotation_colors = annot_colors_ls,
         cellheight = 9,
         cellwidth = 9,
         filename = 'plots/secretome_heatmap_t.pdf'
) 





save.image('rdata/secretome.Rdata')


load('rdata/secretome.Rdata')




significant_cols


cazy_desc = read_tsv('annotation/CAZyDB.07312018.fam-activities.txt', skip = 2, col_names = c('id', 'desc')) %>%
  mutate(type = 'CAZYME')

viru_desc <- read_delim('annotation/virulence_descriptions.txt', col_names = c('id', 'desc'), delim = ',') %>%
  mutate(type = 'VIRULENCE')

card_desc = read_tsv('annotation/aro.tsv') %>%
  mutate(Accession = sub('ARO:', '', Accession)) %>%
  rename('id' = Accession) %>%
  rename('desc' = Name) %>%
  select(-ID)

acces_annot <- read_tsv('annotation/secretome_accessory_annots.txt') %>% select(-X5)

desc_all <- bind_rows(cazy_desc, viru_desc, card_desc)%>% 
  left_join(acces_annot) %>%
  mutate(desc = ifelse(is.na(desc), yes = MainDescription, no = desc)) %>%
  mutate(Description = ifelse(is.na(Description), yes = AccessoryDescription, no = Description)) %>%
  select(-MainDescription, -AccessoryDescription)
 







result_col$hclust$labels[result_col$hclust$order] %>% data.frame(id = .) %>%
  mutate(significant = ifelse(id %in% significant_cols, yes = TRUE, no = FALSE)) %>%
  mutate(heatmap_order = seq_along(significant)) %>%
  left_join(., desc_all, by = 'id') %>%
  mutate(type = ifelse(is.na(type), no = type, yes = case_when(
    grepl('MER', id) ~ 'MEROPS',
    grepl('abH', id) ~ 'LED',
    grepl('300', id) ~ 'CARD',
    TRUE ~ type
  ))) %>% 
  left_join(., gene_representation_df, by = 'id')  %>% 
  write_csv(., path = 'tables/secretome_heatmap_genes.csv')



