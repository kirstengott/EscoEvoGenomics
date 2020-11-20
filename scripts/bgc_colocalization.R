library(GenomicRanges)
library(tidyverse)
library(gggenes)

make_BGC_plot <- function(gff3_file, extra_track = FALSE){
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
  gff <- rtracklayer::import(gff3_file, format = 'gff3') %>%
    data.frame() %>%
    mutate(direction = case_when(
      strand == "+" ~ 1,
      strand == "-" ~ -1,
      TRUE          ~ NA_real_
    )) %>%
    mutate(seqnames = sub(" .*$", "", seqnames)) %>%
    group_by(seqnames, start, end) %>%
    mutate(cluster = sub("genomeID:", "", Note[[1]][which(grepl('genomeID:', Note[[1]]))])) %>%
    ungroup() %>%
    select(-Note, -Parent, -Dbxref) %>%
    mutate(genome = sub("_.*$", "", cluster))
  
  extra_track <- gff %>% filter(type == 'cluster') %>% 
    select(cluster, genome, seqnames) %>%
    inner_join(extra_track, ., by = c('genome', 'seqnames'))
  
  #gff <- bind_rows(gff, extra_track)
  
  product <- sub("product:", "", gff$Note[[1]][which(grepl('prodhuct:', gff$Note[[1]]))])
  
  
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
  
  
    plot_out <- gff %>%
      filter(type %in% c('CDS')) %>%
      ggplot(aes(xmin = start,
                 xmax = end,
                 y = cluster,
                 forward = direction)) +
      geom_gene_arrow(fill = 'white',
                      arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm"),
                      arrow_body_height = unit(3, 'mm'),
                      size = 0.9) +
      geom_gene_arrow(data = extra_track,
                      aes(xmin = start, xmax = end, y = plot_y,
                          forward = direction),
                      fill = 'purple',
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

  if(length(which(grepl('BGCid:', gff$Note[[1]]))) > 0){
    known_bgc <- sub("BGCid:", '', gff$Note[[1]][which(grepl('BGCid:', gff$Note[[1]]))])
    plot_out + labs(subtitle = paste0('Known BGC:', known_bgc))
  } else {
    plot_out
  }
  
}
  
get_card_overlaps <- function(x, annot){
  files = list.files(x, full.names = TRUE)
  card_overlaps <- lapply(files, function(f){
    print(f)
    base = basename(f)
    gff_prefix = strsplit(base, split = '-')[[1]][1]
    
    gff3_file = list.files('antismash_cluster_analysis/', pattern = gff_prefix, full.names = TRUE)
    genome =rev(strsplit(f, split = '-')[[1]])[1]
    terp1 <- read_tsv(f, 
                      col_names = FALSE) 
    terp1 %>%
      select(-X1:-X9) %>%
      filter(X12 == 'CDS') %>%
      mutate(gene = sub("ID=", "", sub(":cds;.*$", "", X18))) %>%
      mutate(genome = genome) %>%
      left_join(., annot, by = c('genome', 'gene')) %>%
      mutate(direction = case_when(
        X16 == "+" ~ 1,
        X16 == "-" ~ -1,
        TRUE          ~ NA_real_
      )) %>%
      filter(## removing things no longer in swissprot
        !is.na(annot)
      )  %>%
      rename('seqnames' = X10,
             'start' = X13,
             'end' = X14,
             'strand' = X16,
             'source' = X11,
             'score' = X15,
             'phase' = X17,
             'Name' = tool) %>%
      select(-score) %>%
      mutate_at(vars(phase), as.numeric) %>%
      mutate(file = base) %>%
      mutate(plot_y = paste0(genome, ":", annot))
  }) %>%
    bind_rows()
}

annot <- read_tsv('annotation/all_annotations.txt', col_names = c('genome', 'gene', 'tool', 'annot')) %>% 
  filter(tool %in% 'CARD')


overlaps_dirs <- dir('antismash_cluster_analysis', pattern = 'output', full.names = TRUE)
overlaps_dirs <- overlaps_dirs[!grepl('pdf', overlaps_dirs)]

x <- overlaps_dirs[1]
y <- cluster_files[1]

lapply(overlaps_dirs, function(x){
  card_overlaps <- get_card_overlaps(x, annot = annot)
  cluster_files <- list.files('antismash_cluster_analysis', pattern = 'gff', full.names = TRUE)
  print(x)
  plots <- lapply(cluster_files, function(y){
    gff3_file <- y
    extra_track = card_overlaps
    outfile_name = paste0(sub(".gff3", "", y), "-", basename(x), ".pdf")
    plot <- make_BGC_plot(gff3_file, extra_track = extra_track)
    ggsave(plot, filename = outfile_name, width = 14, height = 10)
  })
})

