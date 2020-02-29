library(tidyverse)
library(topGO)
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

id_map <- sapply(unique(annots_all$genome), grep, x = unique(all_ortho$genome), value = TRUE)

id_map1 <- names(id_map)
names(id_map1) <- id_map

all_ortho <- all_ortho %>%
  rowwise() %>%
  mutate(genome = id_map1[[genome]]) %>%
  group_by(ortho) %>%
  mutate(ffa_count = length(which(unique(genome) %in% ffa))) %>%
  mutate(trich_count = length(which(unique(genome) %in% trich))) %>%
  mutate(clad_count = length(which(unique(genome) %in% clad)))

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
  select(ortho, annot) %>%
  mutate(all_go = paste(annot, collapse = ",")) %>%
  select(-annot) %>%
  distinct() %>%
  mutate(all_go = paste(unique(strsplit(all_go, ",")[[1]]), collapse = ","))


ortho2go        <- ortho2go_df$all_go
names(ortho2go) <- ortho2go_df$ortho

ortho2go <- lapply(ortho2go, function(x){strsplit(x, ',')[[1]]})

geneNames <- names(ortho2go)
head(geneNames)

myInterestingGenes <- not_ffa
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)

enrich <- lapply(c('BP', 'MF'), function(x){
  GOdata <- new("topGOdata", ontology = x, allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = ortho2go)
  test.stat <- new("classicCount", 
                   testStatistic = GOFisherTest, 
                   name = "Fisher test")
  resultFisher <- getSigGroups(GOdata, test.stat)
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

all_enrich_df <- bind_rows(enrich[[1]][1], enrich[[2]][1])

#printGraph(GOdata, resultFisher, firstSigNodes = nrow(df2), fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
#View(df2)

## Bimap interface:
# Convert the object to a list
an<- c(as.list(GOBPANCESTOR), as.list(GOMFANCESTOR))

off <- c(as.list(GOBPOFFSPRING), as.list(GOMFOFFSPRING))
# Remove GO IDs that do not have any ancestor
an <- an[!is.na(an)]


head(all_enrich_df$go_id)

my_an = an[all_enrich_df$go_id]
#my_an <- my_an[!is.null(my_an)]
my_off = off[all_enrich_df$go_id]

l <- lapply(head(all_enrich_df$go_id, 5), function(x){
  c(all_enrich_df$go_id[grep(x, my_an)], 
    x, 
    all_enrich_df$go_id[grep(x, my_off)])
})

names(l) <- all_enrich_df$go_id

union(l)

reduced_in <- ortho2go_df %>% filter(ortho %in% myInterestingGenes) %>%
  mutate(all_go = strsplit(all_go, ',')) %>%
  unnest(all_go) %>% 
  dplyr::rename('go_id' = all_go) %>%
  right_join(., all_enrich_df, by = 'go_id') %>% 
  group_by(go_id) %>%
  mutate(num_ortho = length(unique(ortho))) %>% 
  dplyr::select(-ortho) %>%
  distinct() %>% 
  filter(num_ortho > 1)

reduced_in$Term = factor(reduced_in$Term, 
                         levels = arrange(reduced_in, num_ortho)$Term)

reduced_in %>%
  ggplot(aes(x = Term, y = num_ortho)) +
  geom_col() +
  coord_flip() +
  theme_minimal() + 
  xlab('GO Term') + 
  ylab('Number of OrthoGroups') + 
  geom_text(aes(y = num_ortho + 35, x = Term, label = signif(pval, 3)), size = 2.5) +
  scale_y_continuous(expand = c(0.12,0)) + 
  theme(axis.title.y = element_text(margin=margin(0,-10,0,0)),
        axis.text.y = element_text(margin = margin(0,-20,0,0)),
        legend.title = element_blank()) +
  ggsave('plots/reduction_GO_enrich.pdf', width = 10)






