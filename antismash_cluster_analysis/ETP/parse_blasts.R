library(tidyverse)


## run in the shell
## for i in `ls ../../gene_annotation/*protein* | grep -v phr | grep -v pin | grep -v psq`; do base=`basename $i`; blastp -query chaetocin.fa -db $i -outfmt 6 -evalue 0.001 -out protein_output/$base; done

blasts <- list.files('antismash_cluster_analysis/ETP/protein_output/', full.names = TRUE)

blast_in <- lapply(blasts, function(x){
  read_tsv(x, col_names = c("query", "subject", "perc_id", 
                                             "align_len", "mismatches", "gaps", 
                                             "q_start", "q_end", "s_start", "s_end", 
                                             "e_val", "bit_score")) %>%
    mutate(genome = basename(x))
})




## 35 percent identity cutoff looks to get rid of most erroneous matches
blast_in %>%
  bind_rows %>%
  group_by(genome, query) %>%
  filter(perc_id >= 35, align_len >= 400) %>% 
  mutate(rank = rank(-bit_score, ties.method="first")) %>%
  filter(rank %in% c(1, 2)) %>% 
  ungroup() %>%
  group_by(genome, subject) %>%
  summarize(mean_percid = mean(perc_id), mean_aln_len = mean(align_len)) %>%
  write_tsv(., path = 'antismash_cluster_analysis/ETP/parsed_blasts.tsv')
  



lapply(blast_in, function(x){
  genome <- unique(x$genome)
  x %>%
    bind_rows %>%
    group_by(genome, query) %>%
    filter(perc_id >= 35, align_len >= 400) %>% 
    mutate(rank = rank(-bit_score, ties.method="first")) %>%
    filter(rank %in% c(1, 2)) %>% 
    select(-genome) %>%
    write_tsv(., paste0('antismash_cluster_analysis/ETP/parsed_output/', genome), col_names = FALSE)
})
