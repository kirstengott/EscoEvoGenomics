library(tidyverse)


source('~/scripts/make_BGC_plot.R')

metadata <- read_csv('metadata_cytoscape.csv')
meta_levels <- read_csv('metadata.csv') %>% .$acc
munge_gca <- function(x){
  paste(strsplit(x, split = "_")[[1]][c(1,2)], collapse = "_")
}

get_cluster <- function(x){
  spliter = strsplit(x, split = "/")[[1]][c(4, 5)]
  cluster_list <- strsplit(spliter[2], split = "\\.")[[1]]
  cluster = cluster_list[which(grepl('cluster', cluster_list))]
  paste0(spliter[1], ":", cluster)
}

files <- list.files("parsed_networks", pattern = 'network', recursive = TRUE, full.names = TRUE)

## pull in the parsed component network files from bigscape output
all_files <- lapply(files, function(x){
  table <- read_tsv(x, col_names = c('component', 'bigscape')) %>%
    left_join(., metadata, by = 'bigscape') %>%
    group_by(component) %>%
    mutate(num_agricultures = length(unique(Agriculture))) %>%
    ungroup() %>%
    filter(num_agricultures > 1)
  #write_csv(table, path = paste0(x, '_parsed.csv'))
  table
})

names(all_files) <- files

#x <- names(all_files)[1]
#y <- unique(table$component)[1]

lapply(names(all_files), function(x){
  table <- all_files[[x]]
  lapply(unique(table$component), function(y){
    print(paste0(x, ':', y))

    sub <- table %>% filter(., component == y) %>%
      group_by(bigscape) %>%
      mutate(antismash = sub(pattern = paste0(acc, "_"), replacement = "", x = bigscape)) %>%
      mutate(antismash_filename = paste0('antismash4/gff3/', acc, '_', antismash, '.gff3')) %>%
      mutate(antismash_filename = sub('ASM430301v1_genomic_', '', antismash_filename)) %>%
      ungroup()
    known_BGC <- sub[which(grepl('BGC', sub$bigscape)), ] %>% .$bigscape


    if(length(known_BGC) > 0){
      sub <- sub[-which(grepl('BGC', sub$bigscape)),]
      t_title = paste('Known BGC: ', known_BGC)
    } else { t_title = ''}
    fi_a <- sub$antismash_filename
    outfile = paste0('bigscape_gff3/', unique(sub$component), "_", sub(" ", "", Sys.time()), ".gff3")
    file.remove(outfile)

    lapply(fi_a, function(x){file.append(file2 = x, file1 = outfile)})

    # print(fi_a)
    # lapply(fi_a, function(x){
    #   dir <- dirname(x)
    #   base = basename(dir)
    #   gene_file <- paste0(dir, '/geneclusters.txt')
    #   command = paste0('~/scripts/antismash2gff3.py ', x, " ", gene_file," ",  base, " >>", outfile)
    #   #write(command, file = 'commands.sh', append = TRUE)
    # })


    plot_out <- paste0('BGC/', basename(x), "_", y, ".pdf")
    cluster_used = fi_a %>% map(., get_cluster) %>% unlist()
    l_out = ceiling(ifelse(length(cluster_used) <= 10, yes = 3, no = length(cluster_used)/5))
    chuncks <- seq(from = 1, to = length(cluster_used), length.out = l_out)
    string_out = ''
    i = 1
    while (i <= length(chuncks) - 1) {
      n = chuncks[i]
      y = chuncks[i+1]
      if(string_out == ""){
        string_out = paste(cluster_used[n:y], collapse = ", ")
      }else{
        n = chuncks[i] + 1
        string_out = c(string_out, paste(cluster_used[n:y], collapse = ", "))
      }
      i = i +1
    }
    sub_label <- paste0('Clusters Used: ',  paste(string_out, collapse = "\n"))
    make_BGC_plot(gff3_file = outfile, genome_levels = meta_levels) +
      labs(title = t_title, subtitle = sub_label) +
      ggsave(plot_out, width = 40, height = 40, units = 'cm')
  })
})


files <- list.files("parsed_networks", pattern = 'network', recursive = TRUE, full.names = TRUE)

all_files <- lapply(files, function(x){
  table <- read_tsv(x, col_names = c('component', 'bigscape')) %>%
    left_join(., metadata, by = 'bigscape') %>%
    group_by(component) %>%
    mutate(num_agricultures = length(unique(Agriculture))) %>%
    ungroup()
  table
})

names(all_files) <- files
all_data <- lapply(files, function(x){
  all_files[[x]] %>%
    mutate(BGC_type = dirname(x))
}) %>% bind_rows() %>%
  mutate(Presence = 1) %>%
  mutate(Agriculture = ifelse(is.na(Agriculture),
                              yes = 'Outgroup',
                              no = Agriculture)) %>%
  filter(!is.na(acc),
         num_agricultures > 0,
         acc != 'g1.txt') %>%
  select(component, BGC_type, acc, Presence)


tree_order <- unique(metadata$tree_order)






# lapply(unique(all_data$BGC_type), function(x){
#   ## per genome, presence absence across all
#   print(x)
#   nrps <- all_data %>% filter(BGC_type == x) %>%
#     distinct() %>%
#     spread(component, Presence) %>%
#     select(-BGC_type) %>%
#     data.frame()
#
#   rownames(nrps) <- nrps$acc
#   nrps$acc <- NULL
#   nrps <- as.matrix(nrps)
#   nrps[is.na(nrps)] <- 0
#
#   rows_order <- tree_order[which(tree_order %in% rownames(nrps))]
#
#   m <- metadata %>% select(Agriculture, acc) %>% distinct() %>% data.frame()
#   rownames(m) <- m$acc
#   m$acc <- NULL
#
#   annotation_row = data.frame(
#     "Agriculture" = m[rownames(nrps), ],
#     row.names = rownames(nrps)
#   )
#   ann_colors = list(
#     'Agriculture' = c(Leafcutter = 'green4', Lower = 'yellow', Coral = 'magenta3', Higher = 'blue3', 'NA' = 'white')
#   )
#   pheatmap::pheatmap(nrps[rows_order,],
#                      #        breaks = breaks,
#                      legend = FALSE,
#                      color = c('grey87', 'black'),
#                      cluster_rows = FALSE,
#                      cluster_cols = TRUE,
#                      border_color = "grey70",
#                      cellwidth = 10,
#                      cellheight = 10,
#                      annotation_row = annotation_row,
#                      annotation_colors = ann_colors,
#                      main = x,
#                      filename = paste0('composite_', x, '.pdf'))
# })



nrps <- all_data %>%
  distinct() %>%
  unite('comp', c(BGC_type, component)) %>%
  spread(comp, Presence) %>%
  data.frame()

rownames(nrps) <- nrps$acc
nrps$acc <- NULL
nrps <- as.matrix(nrps)
nrps[is.na(nrps)] <- 0

rows_order <- tree_order[which(tree_order %in% rownames(nrps))]

m <- metadata %>% select(Agriculture, acc) %>% distinct() %>% data.frame()
rownames(m) <- m$acc
m$acc <- NULL

ann_colors = list(
  'BGC_Type' = c(NRPS = '#1B9E77', PKS1 = '#D95F02', Terpene = '#7570B3', Others = "#E7298A", PKS.NRP = "#66A61E"),
  'Agriculture' = c(Leafcutter = 'green4', Lower = 'yellow', Coral = 'magenta3', Higher = 'blue3', 'NA' = 'white')
)


annotation_col = data.frame(
  "BGC_Type" = sub("\\_.*$", "", colnames(nrps)),
  row.names = colnames(nrps)
)

pheatmap::pheatmap(nrps[rows_order,],
                   legend = FALSE,
                   color = c('grey87', 'black'),
                   cluster_rows = FALSE,
                   cluster_cols = TRUE,
                   border_color = "grey70",
                   cellwidth = 10,
                   cellheight = 10,
                   annotation_row = annotation_row,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   main = "Components present in at least two agricultures",
                   filename = paste0('all_BGC_components_.pdf'))

