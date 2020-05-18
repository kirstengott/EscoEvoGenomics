library(topGO)
library(tidyverse)

source('scripts/color_palettes.R')

xx <- as.list(GOTERM)

metadata <- read_csv('tables/metadata.csv') %>% 
  filter(!is.na(genome_id))


ffa <- metadata %>% filter(!grepl('Outgroup', Agriculture)) %>% .$acc

trich <- metadata %>% filter(grepl('Outgroup2', Agriculture)) %>% .$acc %>%
  sub("\\..*$", "", .)

clad  <- metadata %>% filter(grepl('Outgroup1', Agriculture)) %>% .$acc %>%
  sub("\\..*$", "", .)


annots_all <- read_tsv('annotation/all_annotations.txt', 
                         col_names = c('genome', 'gene', 'tool' ,'annot')) %>%
  filter(genome != 'genome')






load('rdata/orthologues_long.rdata')

all_ortho <- orthog_long %>% dplyr::rename('gene' = genes) %>%
  filter(genome != "GCA_000225605.1_CmilitarisCM01_v01_protein", !is.na(gene)) %>%
  mutate(genome = sub('\\..*$', '', genome))
rm(orthog_long)

id_map         <- sapply(unique(annots_all$genome), grep, x = unique(all_ortho$genome), value = TRUE)
id_map1        <- names(id_map)
names(id_map1) <- id_map


## make sure all of the genomes are in the annotation (will fail if not)
lapply(unique(all_ortho$genome), function(x){
  print(x)
  id_map1[[x]]
})



all_ortho <- all_ortho %>%
  rowwise() %>%
  mutate(genome = id_map1[[genome]]) %>%
  group_by(Orthogroup) %>%
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

annot_sum$Orthogroup %>% n_distinct()


no_out <- filter(annot_sum, 
                 trich_count == 0, 
                 clad_count == 0) %>% 
  .$Orthogroup %>% unique() 


not_ffa <- filter(annot_sum, ffa_count == 0) %>%  .$Orthogroup %>%
  unique()


ortho2go_df <- annot_sum %>% 
  dplyr::select(Orthogroup, annot) %>%
  group_by(Orthogroup) %>%
  mutate(all_go = paste(unique(annot), collapse = ",")) %>%
  dplyr::select(-annot) %>%
  distinct() #%>%
  #mutate(all_go = paste(unique(strsplit(all_go, ",")[[1]]), collapse = ","))


ortho2go        <- ortho2go_df$all_go
names(ortho2go) <- ortho2go_df$Orthogroup

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


reduced_in <- ortho2go_df %>% filter(Orthogroup %in% myInterestingGenes) %>%
  mutate(all_go = strsplit(all_go, ',')) %>%
  unnest(all_go) %>%
  dplyr::rename('go_id' = all_go) %>%
  right_join(., enriched_leaves, by = 'go_id') %>%
  left_join(., all_ortho, by = 'Orthogroup') %>%
  group_by(go_id, clade_groups) %>%
  mutate(num_ortho = case_when(
    clade_groups == 'Trichoderma' ~ (length(Orthogroup)/nrow(filter(metadata, genus == 'Trichoderma'))),
    clade_groups == 'Hypomyces_Cladobotryum' ~ (length(Orthogroup)/nrow(filter(metadata, clade_groups == 'Hypomyces_Cladobotryum')))
  )) %>%
  mutate(ortho_pass = ifelse(num_ortho >= 10, TRUE, FALSE)) %>%
  distinct() %>%
  filter(any(ortho_pass)) %>% 
  dplyr::select(go_id, clade_groups, Term, Ontol, BH.correction, num_ortho) %>% 
  distinct() %>%
  ungroup() %>%
  filter(Ontol != 'CC') %>%
  mutate(Ontol = case_when(
    Ontol == 'MF' ~ 'Molecular Function',
    Ontol == 'BP' ~ 'Biological Process'
  )) %>%
  mutate(Term_s = case_when(
    Term == "oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen" ~ "oxidoreductase activity",
    Term == "DNA-binding transcription factor activity, RNA polymerase II-specific" ~ "DNA-binding transcription factor activity",
    TRUE ~ Term
  )) %>%
  mutate(pval_s = case_when(
    BH.correction <= 0.0001  ~ "****",
    BH.correction <= 0.001 ~ "***",
    BH.correction <= 0.01 ~ "**",
    BH.correction <= 0.5 ~ "*",
    TRUE ~ "NA"
  ))

pvals <- reduced_in %>% filter(clade_groups == 'Trichoderma') %>%
  dplyr::select(-clade_groups) %>% distinct()

reduced_in$Term_s = factor(reduced_in$Term_s,
                         levels = unique(arrange(reduced_in, num_ortho)$Term_s))

write_csv(reduced_in, path = 'tables/reduced_go_terms.csv')

fig <- reduced_in %>%
  ggplot(aes(x = Term_s, y = num_ortho)) +
  geom_col(position = 'dodge', aes(fill = clade_groups)) +
  theme_classic() +
  xlab('GO Term') +
  coord_flip() +
  ylab('Average Number of OrthoGroups') +
  ggtitle('GO term Enrichment of Orthogroups Lost in Escovopsis') +
  geom_text(data = pvals, aes(y = num_ortho + 20, x = Term_s, label = pval_s, angle = 90), size = 2.5) +
  facet_wrap(~Ontol, nrow = 2, ncol = 1, scales = 'free_y', strip.position = 'left') +
  scale_y_continuous(expand = c(0,0), limits = c(0, max(reduced_in$num_ortho)+40),
                     breaks = seq(0, max(reduced_in$num_ortho)+40, by = 50),
                     labels = seq(0, max(reduced_in$num_ortho)+40, by = 50)
                     ) +
  theme(#axis.title.y = element_text(margin=margin(0,-10,0,0)),
        #axis.text.y = element_text(margin = margin(0,-20,0,0)),
    panel.border = element_rect(colour = '#000000', fill = NA, size = 1),
    legend.title = element_blank(),
    legend.position = "none")   +
  scale_fill_manual(values = clade_colors)
fig

fig + ggsave('plots/reduction_GO_enrich.pdf', width = 7)






