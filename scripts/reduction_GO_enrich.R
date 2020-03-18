library(topGO)
library(tidyverse)

xx <- as.list(GOTERM)

metadata <- read_csv('tables/metadata.csv')


ffa <- metadata %>% filter(!grepl('Outgroup', Agriculture)) %>% .$acc

trich <- metadata %>% filter(grepl('Outgroup2', Agriculture)) %>% .$acc %>%
  sub("\\..*$", "", .)

clad  <- metadata %>% filter(grepl('Outgroup1', Agriculture)) %>% .$acc %>%
  sub("\\..*$", "", .)


annots_all <- read_tsv('annotation/all_annotations.txt', 
                         col_names = c('genome', 'gene', 'tool' ,'annot'))

all_ortho <- read_csv('annotation/fastortho/orthologues_full.csv',
                      col_names = c('ortho', 'genome', 'gene'))
id_map         <- sapply(unique(annots_all$genome), grep, x = unique(all_ortho$genome), value = TRUE)
id_map1        <- names(id_map)
names(id_map1) <- id_map
all_ortho <- all_ortho %>%
  rowwise() %>%
  mutate(genome = id_map1[[genome]]) %>%
  group_by(ortho) %>%
  mutate(ffa_count = length(which(unique(genome) %in% ffa))) %>%
  mutate(trich_count = length(which(unique(genome) %in% trich))) %>%
  mutate(clad_count = length(which(unique(genome) %in% clad))) %>% 
  group_by(genome) %>%
  mutate(acc = grep(genome, metadata$acc, value = TRUE)) %>%
  ungroup()%>%
  left_join(., metadata, by = 'acc')


joined_ortho <- all_ortho %>%
  left_join(annots_all, by = c('genome', 'gene'))


#dplyr::filter(joined_ortho, tool %in% c('LED', 'MEROPS', 'cazy'), ffa_count == 0) %>% 
#  .$ortho %>% unique() 
  

annot_sum <- joined_ortho %>%
  filter(tool == 'GO') 

annot_sum$ortho %>% n_distinct()


no_out <- filter(annot_sum, 
                 trich_count == 0, 
                 clad_count == 0) %>% 
  .$ortho %>% unique() 


not_ffa <- filter(annot_sum, ffa_count == 0) %>%  .$ortho %>%
  unique()


ortho2go_df <- annot_sum %>% 
  dplyr::select(ortho, annot) %>%
  mutate(all_go = paste(annot, collapse = ",")) %>%
  dplyr::select(-annot) %>%
  distinct() %>%
  mutate(all_go = paste(unique(strsplit(all_go, ",")[[1]]), collapse = ","))


ortho2go        <- ortho2go_df$all_go
names(ortho2go) <- ortho2go_df$ortho

ortho2go <- lapply(ortho2go, function(x){strsplit(x, ',')[[1]]})

geneNames <- names(ortho2go)
head(geneNames)



## define the genes that I'm interested in (those not in fungus farming ant escovopsis)
myInterestingGenes <- not_ffa
geneList           <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList)    <- geneNames


## perform an enrichment test on the genes that are not in FFA
enrich <- lapply(c('BP', 'MF', 'CC'), function(x){
  GOdata <- new("topGOdata", ontology = x, allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = ortho2go)
  test.stat <- new("classicCount", 
                   testStatistic = GOFisherTest, 
                   name = "Fisher test")
  resultFisher         <- getSigGroups(GOdata, test.stat)
  pVal                 <- data.frame(pval=signif(score(resultFisher), 6),BH.correction=seq(1, length(score(resultFisher))))
  pVal                 <- pVal[order(pVal$pval),]
  pVal$BH.correction   <- signif(p.adjust(pVal$pval, method="BH"), 6)
  pVal.sub             <- pVal[pVal$BH.correction<0.01,]
  pVal.sub$go_id       <- rownames(pVal.sub)
  pVal.try <- cbind(pVal.sub, Term=sapply(pVal.sub$go_id, FUN=function(n){Term(xx[[n]])}), Ontol=sapply(pVal.sub$go_id, FUN=function(n){Ontology(xx[[n]])}))
  
  if (nrow(pVal.try) >0){
    df2                  <- pVal.try[,c("go_id", "Term", "Ontol", "pval", "BH.correction")]
    out <- list(df2, GO, resultFisher)
    invisible(out)
  } else {
    message("No enrichment")
  }
})



all_enrich_df <- bind_rows(enrich[[1]][1], enrich[[2]][1], enrich[[3]][1])
all_enrich_df %>% dplyr::select(go_id, BH.correction, everything()) %>% write_tsv('tables/ffa_reduction_GO.tsv')


#printGraph(GOdata, resultFisher, firstSigNodes = nrow(df2), fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#View(df2)

## get lists of ancestors and offspring
an  <- c(as.list(GOBPPARENTS), as.list(GOMFPARENTS), as.list(GOCCPARENTS))
off <- c(as.list(GOBPCHILDREN), as.list(GOMFCHILDREN), as.list(GOCCCHILDREN))
# Remove GO IDs that do not have any ancestor
an <- an[!is.na(an)]

my_an <- unlist(an[all_enrich_df$go_id])

enriched_leaves <- all_enrich_df %>%
  filter(!go_id %in% my_an)


reduced_in <- ortho2go_df %>% filter(ortho %in% myInterestingGenes) %>%
  mutate(all_go = strsplit(all_go, ',')) %>%
  unnest(all_go) %>%
  dplyr::rename('go_id' = all_go) %>%
  right_join(., enriched_leaves, by = 'go_id') %>%
  left_join(., all_ortho, by = 'ortho') %>%
  dplyr::rename('OrthoGroup' = ortho) %>%
  group_by(go_id, clade_groups) %>%
  mutate(num_ortho = case_when(
    clade_groups == 'Trichoderma' ~ (length(OrthoGroup)/nrow(filter(metadata, genus == 'Trichoderma'))),
    clade_groups == 'Hypomyces/Cladobotryum' ~ (length(OrthoGroup)/nrow(filter(metadata, clade_groups == 'Hypomyces/Cladobotryum')))
  )) %>%
  mutate(ortho_pass = ifelse(num_ortho >= 10, TRUE, FALSE)) %>%
  distinct() %>%
  filter(any(ortho_pass)) %>% 
  select(go_id, clade_groups, Term, Ontol, BH.correction, num_ortho) %>% 
  distinct() %>%
  ungroup()

pvals <- reduced_in %>% filter(clade_groups == 'Trichoderma') %>%
  select(-clade_groups) %>% distinct()

reduced_in$Term = factor(reduced_in$Term,
                         levels = unique(arrange(reduced_in, num_ortho)$Term))
reduced_in %>%
  ggplot(aes(x = Term, y = num_ortho)) +
  geom_col(position = 'dodge', aes(fill = clade_groups)) +
  coord_flip() +
  theme_minimal() +
  xlab('GO Term') +
  ylab('Average Number of OrthoGroups') +
  ggtitle('GO term Enrichment of Orthogroups Lost in Escovopsis') +
  geom_text(data = pvals, aes(y = num_ortho + 35, x = Term, label = signif(BH.correction, 3)), size = 2.5) +
  scale_y_continuous(expand = c(0.12,0)) +
  facet_wrap(~Ontol, nrow = 3, ncol = 1, scales = 'free_y') +
  theme(axis.title.y = element_text(margin=margin(0,-10,0,0)),
        axis.text.y = element_text(margin = margin(0,-20,0,0)),
        legend.title = element_blank())   +
  scale_fill_manual(values = c('Trichoderma' = 'black', 'Hypomyces/Cladobotryum' = 'gray'))+
  ggsave('plots/reduction_GO_enrich.pdf', width = 10)






