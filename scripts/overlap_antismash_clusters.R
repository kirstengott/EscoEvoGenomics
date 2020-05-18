
library(GenomicRanges)
library(tidyverse)





parsed_networks <- list.files('bigscape/parsed_networks/', full.names = TRUE) ## don't look at the mix

net_in <- lapply(parsed_networks, function(x){
  read_tsv(x, col_names = c('component', 'id')) %>%
             mutate(type = sub('_c0.30.network', '', basename(x)))
}) %>% bind_rows() %>%
  mutate(source = case_when(
    grepl('maker-glimmerHMM', id) ~ 'maker;glimmerHMM',
    grepl('maker', id) ~ 'maker',
    grepl('glimmerHMM', id) ~ 'glimmerHMM',
    TRUE ~ 'NA'
  )) %>%
  filter('source' != 'NA') %>%
  group_by(id)%>%
  mutate(cluster_id = rev(strsplit(id, split = '\\.')[[1]])[1])%>%
  mutate(id_keep = paste0(type, ".", component)) %>%
  mutate(genome = sub("\\..*$", "", id)) %>%
  mutate(genome = sub("_maker.*$", "", genome)) %>%
  mutate(genome = sub("_glimmerHMM.*$", "", genome)) %>%
  ungroup() %>%
  select(source, cluster_id, id_keep, genome) %>%
  filter(!grepl('BGC', genome))

   
as_clusters <- read_csv('tables/antismash_clusters_supplement.csv', col_names = c("chr", 
                                                                                  'start', 
                                                                                  'stop', 
                                                                                  'strand', 
                                                                                  'source', 
                                                                                  'maker_id', 
                                                                                  'glimmer_id', 
                                                                                  'acc')) %>%
mutate(cluster_id = ifelse(is.na(maker_id), yes = glimmer_id, no = maker_id)) %>%
  mutate(genome = sub("\\..*$", '', acc)) %>%
  left_join(., net_in, by = c('genome', 'source', 'cluster_id'))


gff1        <- GRanges(as_clusters)





gene_map <- read_delim('tables/trichoderma_gene_map.txt', delim = ' ', col_names = c('gene_id', 'prot_id'))

gff <- list.files('gene_annotation', full.names = TRUE)

mine <- c("gene_annotation/GCA_004303015.1_ASM430301v1_genomic.genemodels.gff",
          "gene_annotation/GCA_011799845.1.genemodels.gff",
          "gene_annotation/ICBG1054.genemodels.gff",
          "gene_annotation/ICBG1065.genemodels.gff",
          "gene_annotation/ICBG1075.genemodels.gff",
          "gene_annotation/ICBG1096.genemodels.gff",
          "gene_annotation/ICBG2046.genemodels.gff",
          "gene_annotation/ICBG2047.genemodels.gff",
          "gene_annotation/ICBG2048.genemodels.gff",
          "gene_annotation/ICBG2049.genemodels.gff",
          "gene_annotation/ICBG710.genemodels.gff",
          "gene_annotation/ICBG712.genemodels.gff",
          "gene_annotation/ICBG721.genemodels.gff",
          "gene_annotation/ICBG726.genemodels.gff",
          "gene_annotation/ICBG730.genemodels.gff",
          "gene_annotation/ICBG731.genemodels.gff",
          "gene_annotation/ICBG733.genemodels.gff",
          "gene_annotation/ICBG736.genemodels.gff",
          "gene_annotation/ICBG742.genemodels.gff",
          "gene_annotation/ICBG751.genemodels.gff",
          "gene_annotation/NIGB00000000.genemodels.gff",
          "gene_annotation/NIGC00000000.genemodels.gff",
          "gene_annotation/NIGD00000000.genemodels.gff",
          "gene_annotation/NQYQ00000000.genemodels.gff",
          "gene_annotation/NQYR00000000.genemodels.gff",
          "gene_annotation/NQYS00000000.genemodels.gff",
          "gene_annotation/SPDT00000000.1_genomic.genemodels.gff")
x <- "gene_annotation/SPDT00000000.1_genomic.genemodels.gff"
gff2 <- lapply(gff, function(x){
  print(x)
  t <- rtracklayer::import(x, format = 'gff3')
  m <- mergeByOverlaps(gff1, t) %>%
    data.frame(.) %>% filter(type == 'mRNA')
  if(x %in% mine){
    out <- m %>% select(Name, id_keep)  %>% distinct() %>%
      rename('prot_id' = Name)
  } else if(x == 'gene_annotation/LGSR00000000.genemodels.gff'){
    gid <- m %>% select(orig_protein_id, id_keep) %>% mutate(orig_protein_id = sub("-.*$", "", strsplit(orig_protein_id, split = '\\|')[[1]][3])) 
    out <- gid %>% rename('gene_id' = orig_protein_id) %>%
                  left_join(., gene_map,by = 'gene_id') %>% select(-gene_id) %>% distinct()
  } else {
    out <- m %>% select(locus_tag, id_keep) %>% distinct() %>%
      rename('gene_id' = locus_tag) %>%
      left_join(., gene_map, by = 'gene_id') %>%
      select(-gene_id) %>%
      distinct()
  }
  out %>%
    mutate(genome = sub("\\..*$", "", basename(x)))
  }) %>%
  bind_rows()


write_csv(gff2, path = 'tables/antismash_genes.csv')









