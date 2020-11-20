library(GenomicRanges)
library(tidyverse)
library(gggenes)

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

color_mapping <- c('ACP' =  "#1a476f",
                   'AMP-binding' = "#90353b",
                   'Aminotran_1_2' = "#55752f",
                   'Condensation'= "#e37e00",
                   'ECH' = "#6e8e84",
                   'MT' = "#c10534",
                   'NAD_binding_4' = "#938dd2",
                   'PCP' = "#cac27e",
                   'PKS_AT' = "#a0522d",
                   'PKS_DH' = "#7b92a8",
                   'ER' = "#2d6d66",
                   'PKS_KR' = "#9c8847",
                   'PKS_KS' = "#bfa19c",
                   'TD' = "#ffd200",
                   'Thioesterase' = "#d9e6eb",
                   'Additional Biosynthetic Domains' = "#404040",
                   'virulence' = 'purple',
                   'CARD' = 'blue')


#x <- files[1]

gff3_file <- 'antismash_cluster_analysis/ETP/gff_full/all.gff3'


gff <- rtracklayer::import(gff3_file, format = 'gff3') %>%
  data.frame() %>%
  mutate(direction = case_when(
    strand == "+" ~ 1,
    strand == "-" ~ -1,
    TRUE          ~ NA_real_
  )) %>%
  mutate(seqnames = sub(" .*$", "", seqnames))


clusters <- lapply(seq(1, nrow(gff)), function(x){
  ls = unlist(gff[x, 'Note'])
  ls[which(grepl('GenomeID:', ls, ignore.case = TRUE))] %>% 
    sub("GenomeID:", "", .) %>%
    sub("genomeID:", "", .)
  }) %>% unlist()


gff$clusters <- clusters

genomes <- lapply(clusters, function(x){
  x %>% sub("_maker.*$", "", .) %>% sub("\\..*$", "", .)
}) %>% unlist()

gff$genome <- genomes

#product <- sub("product:", "", gff$Note[[1]][which(grepl('prodhuct:', gff$Note[[1]]))])


## may integrate later somehow, not sure how
core_biosynthetic <- gff %>%
  filter(type == 'aSDomain') %>%
  mutate(ID = sub("PKS_", "", ID)) %>%
  GRanges()

## may need to make more consice later
subdomains <- gff %>% 
  filter(type == 'PFAM_domain') %>%
  GRanges()

overlaps <- mergeByOverlaps(query = subdomains, subject = core_biosynthetic) %>%
  data.frame()


sd_data <-  gff %>% filter(type == 'PFAM_domain') %>%
  mutate(color = ifelse(Name %in% overlaps$subdomains.Name,
                        yes = overlaps$core_biosynthetic.Name,
                        no = 'Additional Biosynthetic Domains')) %>%
  mutate(color = ifelse(Name %in% unique(core_biosynthetic$Name),
                        yes = Name,
                        no = color)) #%>%
# mutate(color = ifelse(Name == 'virulence', yes = 'virulence', no = color)) %>%
# mutate(color = ifelse(Name == 'CARD', yes = 'CARD', no = color))

gff$genome <- factor(gff$genome, levels = order)

plot_out <- gff %>%
  filter(type %in% c('CDS'))  %>%
  ggplot(aes(xmin = start,
             xmax = end,
             y = clusters,
             forward = direction)) +
  geom_gene_arrow(fill = 'white',
                  arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(1, "mm"),
                  arrow_body_height = unit(3, 'mm'),
                  size = 0.9) +
  facet_wrap(~ genome, scales = "free_y", ncol = 1) +
  geom_gene_arrow(data = sd_data,
                  arrowhead_height = unit(0, "mm"),
                  arrowhead_width = unit(0, "mm"),
                  arrow_body_height = unit(2.1, 'mm'),
                  aes(fill = color)) +
  scale_fill_manual(values = color_mapping) +
  geom_gene_arrow(fill = NA,
                  arrowhead_height = unit(3, "mm"),
                  arrowhead_width = unit(1, "mm"),
                  arrow_body_height = unit(3, 'mm'),
                  size = 0.9) +
  theme_genes()


ggsave(plot_out, filename = 'ETP.pdf', width = 14, height = 15)




gff %>% filter(type == 'cluster') %>% select(genome, clusters) %>%
  distinct() %>%
  View()

# if(length(which(grepl('BGCid:', gff$Note[[1]]))) > 0){
#   known_bgc <- sub("BGCid:", '', gff$Note[[1]][which(grepl('BGCid:', gff$Note[[1]]))])
#   plot_out + labs(subtitle = paste0('Known BGC:', known_bgc))
# } else {
#   plot_out
# }


#load('antismash/gff_rimg/GCA_000167675.2')
#etp  <- rtracklayer::import('antismash_cluster_analysis/ETP/all_blast_gff.gff')

#all <- mergeByOverlaps(gff, etp)


# overlaps_dirs <- dir('antismash_cluster_analysis/ETP/gff_look', pattern = 'mix', full.names = TRUE)
# 
# 
# lapply(overlaps_dirs, function(x){
#     plot <- make_BGC_plot(gff3_file = x)
#     outfile_name = sub(
#       pattern = "gff3$", replacement = "pdf", x
#     )
#     ggsave(plot, filename = outfile_name, width = 14, height = 10)
#   })



