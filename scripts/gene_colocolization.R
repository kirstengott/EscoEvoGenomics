library(tidyverse)

## test if annotated intervals are related spatially

dir.create('bedtools_fischer')

tool_sub <- c('virulence', 'led', 'MEROPS', 'cazy', 'CARD', 'busco', 'antismash')

annot <- read_tsv("annotation/all_annotations.txt", col_names = c('genome', 'gene', 'tool', 'annot')) %>%
  filter(tool %in% tool_sub)


genome_files <- list.files('gene_annotation', pattern = 'mcscan.gff', full.names = TRUE)

#x <- genome_files[28]
genomes = lapply(genome_files, function(x){
  gff <- read_tsv(x, col_names = c('chr', 'start', 'end', 'gene'))
  genome = basename(x)
  if(grepl('GCA', genome)){
    genome1 = paste(strsplit(genome, "_")[[1]][c(1,2)], collapse = "_") %>% sub("\\..*$", "", .)
  } else {
    genome1 = basename(x) %>% sub("_genomic.genemodels.mcscan.gff", "", .)  %>% sub(".mcscan.gff", "", .) %>% sub("\\..*$", "", .)
  }
  if(genome1 %in% annot$genome){
    lapply(tool_sub, function(y){
      annot_sub <- annot %>% filter(genome == genome1, tool == y)
      out_f <- paste0('bedtools_fischer/', genome1, "-", y, ".bed")
      gff %>% filter(gene %in% annot_sub$gene) %>% write_tsv(., path = out_f, col_names = FALSE)
      g_file = list.files('fasta', pattern = 'names', full.names = TRUE) %>% grep(genome1, ., value = TRUE) %>% grep('cds', ., invert = TRUE, value = TRUE)
      command = paste0('/Users/kirstengotting/miniconda3/bin/sortBed -i ', out_f, " -g ", g_file,  " >", out_f, ".sort")
      system(command = command, intern = FALSE)
    })
    print(genome1)
    return(genome1)
  }
})

genomes <- unlist(genomes)



## add in a 10% overlap
## subtract those results from a larger window (1kb)
#x <- genomes[1]
all_sig_same = lapply(genomes, function(x){
  files <- list.files('bedtools_fischer', pattern = x, full.names = TRUE) %>% grep('sort', ., value = TRUE) %>%
    grep('window', ., value = TRUE, invert = TRUE)
  g_file = list.files('fasta', pattern = 'lens', full.names = TRUE) %>% 
    grep(x, ., value = TRUE) %>% 
    grep('cds', ., invert = TRUE, value = TRUE)
  files_len = length(files)
  for (i in files) {
    if(length(files) == 1){
      break
    }
    dfs <- data.frame(pvalues = as.numeric(),
                      f_test = as.character(),
                      comparison = as.character(),
                      genome = as.character(),
                      stringsAsFactors = FALSE)
    dfs_t <- lapply(seq(2, length(files)), function(y) {
      file1 = i
      file2 = files[y]
      test = paste0(sub("^[^-]*-", "", sub(".bed.sort", "", basename(file1))), "-", sub("^[^-]*-", "", sub(".bed.sort", "", basename(file2))))
      command = paste0("/Users/kirstengotting/miniconda3/bin/bedtools fisher -f 0.1 -F 0.1 -a ", file1, " -b ", file2, " -g ", g_file)
      out <- system(command = command, intern = TRUE)
      out_df <- data.frame(pvalues = out[length(out)] %>% gsub("\t", ",", .) %>%  strsplit(., ",") %>% unlist() %>% as.numeric(), 
               f_test = out[length(out)-1] %>% gsub("\t", ",", .) %>%  strsplit(., ",") %>% unlist(),
               comparison = test,
               genome = sub("-.*$", "", basename(file1)),
               stringsAsFactors = FALSE
               ) %>% filter(pvalues <= 0.05, pvalues != '-nan')
    }) %>% bind_rows()
    dfs <- rbind(dfs, dfs_t)
    files = files[-which(files == i)]
    return(dfs)
  }
}) %>% bind_rows()



all_sig_diff = lapply(genomes, function(x){
  files <- list.files('bedtools_fischer', pattern = x, full.names = TRUE) %>% grep('sort', ., value = TRUE) %>%
    grep('window', ., value = TRUE, invert = TRUE)
  g_file = list.files('fasta', pattern = 'lens', full.names = TRUE) %>% 
    grep(x, ., value = TRUE) %>% 
    grep('cds', ., invert = TRUE, value = TRUE)
  files_len = length(files)
  for (i in files) {
    if(length(files) == 1){
      break
    }
    dfs <- data.frame(pvalues = as.numeric(),
                      f_test = as.character(),
                      comparison = as.character(),
                      genome = as.character(),
                      stringsAsFactors = FALSE)
    dfs_t <- lapply(seq(2, length(files)), function(y) {
      file1 = i
      file2 = files[y]
      test = paste0(sub("^[^-]*-", "", sub(".bed.sort", "", basename(file1))), "-", sub("^[^-]*-", "", sub(".bed.sort", "", basename(file2))))
      file1_window = paste0(file1, ".window")
      file2_window = paste0(file2, ".window")
      file1_window_s = paste0(file1, ".window.sort")
      file2_window_s = paste0(file2, ".window.sort")
      
      command = paste0("/Users/kirstengotting/miniconda3/bin/bedtools slop -b 10000 -i ", file1, " -g ", g_file, " >", file1_window)
      system(command = command, intern = FALSE)
      command = paste0("/Users/kirstengotting/miniconda3/bin/bedtools slop -b 10000 -i ", file2, " -g ", g_file, " >", file2_window)
      system(command = command, intern = FALSE)
      
      command = paste0('/Users/kirstengotting/miniconda3/bin/sortBed -i ', file1_window, " -g ", g_file,  " >", file1_window_s)
      system(command = command, intern = FALSE)
      command = paste0('/Users/kirstengotting/miniconda3/bin/sortBed -i ', file2_window, " -g ", g_file,  " >", file2_window_s)
      system(command = command, intern = FALSE)
      
      
      command = paste0("/Users/kirstengotting/miniconda3/bin/bedtools fisher -a ", file1_window_s, " -b ", file2_window_s, " -g ", g_file)
      out <- system(command = command, intern = TRUE)
      
      out_df <- data.frame(pvalues = out[length(out)] %>% gsub("\t", ",", .) %>%  strsplit(., ",") %>% unlist() %>% as.numeric(), 
                           f_test = out[length(out)-1] %>% gsub("\t", ",", .) %>%  strsplit(., ",") %>% unlist(),
                           comparison = test,
                           genome = sub("-.*$", "", basename(file1)),
                           stringsAsFactors = FALSE
      )
    }) %>% bind_rows()
    dfs <- rbind(dfs, dfs_t)
    files = files[-which(files == i)]
    return(dfs)
  }
}) %>% bind_rows()



metadata <- read_csv('tables/metadata.csv') %>%
  mutate(genome = sub("\\..*$", "", acc)) %>%
  select(genome, Agriculture, genus_species, genus)


all_sig <- anti_join(all_sig_diff, all_sig_same, by = c('f_test', 'comparison', 'genome')) %>% filter(pvalues <= 0.05) %>%
  left_join(., metadata, by = 'genome')

table(all_sig$comparison)
table(all_sig$genome)


all_sig %>% write_csv('tables/gene_colocalization.csv')

all_sig <- read_csv('tables/gene_colocalization.csv')


all_sig %>%
  select(comparison, genome, Agriculture) %>% 
  distinct() %>%
  group_by(Agriculture, comparison) %>%
  summarize(sum = length(unique(genome))) %>% 
  arrange(Agriculture, comparison) %>%
  View()



