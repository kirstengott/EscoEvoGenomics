#library(KEGGREST)
#library(GO.db)
#library(GSEABase)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)
library(UpSetR)
library(tidyverse)


source('scripts/color_palettes.R')

annot <- read_tsv('annotation/all_annotations.txt', col_names = c('acc', 'gene', 'tool', 'annot'))

meta <- read_csv('tables/metadata.csv')


all_gene_id_files <- list.files('annotation/gene_ids/', full.names = TRUE)

all_gene_ids <- lapply(all_gene_id_files, function(x){
  ids <- scan(x, what = 'character')
  data.frame(gene = ids, genome = rep(basename(x), length(ids)), 
             stringsAsFactors = FALSE) %>%
    mutate(acc = case_when(
      grepl('GCA_004303015', genome) ~ 'GCA_004303015.1',
      grepl('GCA_011799845', genome) ~ 'GCA_011799845.1',
      grepl('GCA', genome) ~ paste(strsplit(genome, split = "_")[[1]][c(1,2)], collapse = "_"),
      grepl("SPDT00000000.1_genomic", genome) ~ "SPDT00000000",
      grepl('LGSR', genome) ~ "LGSR00000000",
      grepl('.fasta_gene_ids', genome) ~ sub(".all.maker.proteins.fasta_gene_ids", '', genome),
      TRUE ~ genome
    )) 
}) %>% bind_rows() %>%
  dplyr::select(-genome)


orthog      <- read_tsv('tables/Orthogroups.tsv')
load('rdata/orthologues_long.rdata')




ortho_full  <- orthog_long %>% 
  filter(!is.na(genes)) %>%
  group_by(genome) %>%
  mutate(acc = case_when(
      grepl('GCA_004303015', genome) ~ 'GCA_004303015.1',
      grepl('GCA_011799845', genome) ~ 'GCA_011799845.1',
      grepl('GCA', genome) ~ paste(strsplit(genome, split = "_")[[1]][c(1,2)], collapse = "_"),
      grepl("SPDT00000000.1_genomic", genome) ~ "SPDT00000000",
      grepl('all.maker', genome) ~ sub(".all.maker.proteins", '', genome),
      grepl('LGSR', genome) ~ "LGSR00000000",
      TRUE ~ genome
    )) %>%
  ungroup() %>%
  rename('gene' = genes) %>%
  full_join(all_gene_ids, ., by = c('acc', 'gene')) %>%
  mutate(orthology_meta = ifelse(is.na(Orthogroup),
                                 yes = gene,
                                 no = Orthogroup)) %>%
  left_join(., meta, by = 'acc') %>%
  filter(!is.na(clade_groups)) %>%
  group_by(clade_groups, orthology_meta) %>%
  mutate(n_genomes = length(unique(acc))) %>%
  mutate(n_copies_per_genome = length(acc)) %>%
  ungroup()



## compare the three pangenomes of the three clades

groups <- c('Trichoderma', 'Hypomyces_Cladobotryum', 'Escovopsis')

ag <- c("Lower",
         "Coral",
         "Higher",
         "Leafcutter")


venn1 <- lapply(groups, function(x){
  if(x %in% c('Escovopsis', 'Trichoderma')){
    cutoff <- trunc(length(which(meta$clade_groups %in% x)) *.95)

  } else{
    cutoff <- length(which(meta$clade_groups %in% x))
  }
  ortho_full %>% 
    filter(clade_groups %in% x) %>%
    filter(n_genomes >= cutoff)
}) %>% bind_rows() %>%
  select(orthology_meta, clade_groups, n_genomes) %>%
  distinct() %>%
  group_by(clade_groups, orthology_meta) %>%
  summarize(genome_count = sum(n_genomes)) %>%
  ungroup()


venn <- lapply(ag, function(x){
  cutoff <- trunc(length(which(meta$Agriculture %in% x)) *.95)
  ortho_full %>% 
    filter(Agriculture %in% x) %>%
    group_by(orthology_meta) %>%
    mutate(n_genomes = length(unique(genus_species)))%>%
    ungroup() %>%
    filter(n_genomes >= cutoff)
}) %>% bind_rows() %>%
  select(orthology_meta, Agriculture, n_genomes) %>%
  distinct() %>%
  group_by(Agriculture, orthology_meta) %>%
  summarize(genome_count = sum(n_genomes)) %>%
  ungroup() %>%
  rename('clade_groups' = Agriculture) %>%
  bind_rows(., venn1)


write_csv(venn, path = 'tables/clade_orthogroup_venn.csv')



listInput <- lapply(c(groups, ag), 
                    function(x){
                      venn %>% filter(clade_groups == x) %>% .$orthology_meta %>% unique()
                    })

names(listInput) <- c(groups, ag)

listInput[['AllOutgroup']] <- c(listInput$Trichoderma, listInput$Hypomyces_Cladobotryum)

# listInput <- list('Escovopsis' = esco_set,, 
#                   'Hypomyces/Cladobotryum' = ch_set,
#                   'Trichoderma' = t_set)
# set_name(listInput)

lapply(listInput, head)

pdf("plots/orthologues_upset_summary.pdf", width = 8, height = 8)
upset(fromList(listInput), order.by = "freq", 
      sets = c(groups),
       point.size = 3.5, 
       line.size = 1.5, 
       text.scale = c(1.3, 1.3, 1, 1, 2, 0.75),
       mb.ratio = c(0.8, 0.2))
dev.off()

#pdf("plots/orthologues_upset_subset.pdf", width = 8, height = 8)

#list('Lower','Coral','Higher','Leafcutter','AllOutgroup')

pdf("plots/orthologues_upset_subset.pdf", width = 8, height = 8)
upset(fromList(listInput), order.by = "freq", 
      sets = c(ag, 'AllOutgroup'),
      intersections = list(
        list('Lower','Higher','Leafcutter','AllOutgroup'),
        list('Higher', "Leafcutter"),
        list('Lower', 'AllOutgroup'),
        list('Coral', 'Higher', 'Lower', 'AllOutgroup'),
        list('Coral', 'Higher', 'Leafcutter', 'AllOutgroup'),
        list('Coral', 'Lower', 'Leafcutter', 'AllOutgroup'),
        list('Lower','Coral', 'AllOutgroup'),
        list('Higher','Leafcutter','AllOutgroup'),
        list('Leafcutter'),
        list('Lower'),
        list('Lower','Higher','AllOutgroup'),
        list('Lower','Coral','Higher','Leafcutter'),
        list('Lower','Leafcutter','AllOutgroup')
      ),
      point.size = 3.5, 
      line.size = 1.5, 
      text.scale = c(1.3, 1.3, 1, 1, 2, 0.75),
      mb.ratio = c(0.8, 0.2))
dev.off()

pdf("plots/orthologues_upset_all.pdf", width = 8, height = 8)
upset(fromList(listInput), order.by = "freq", 
      sets = c(ag, 'AllOutgroup'),
      point.size = 3.5, 
      line.size = 1.5, 
      text.scale = c(1.3, 1.3, 1, 1, 2, 0.75),
      mb.ratio = c(0.8, 0.2))
dev.off()
 

 ## Keep this as it show a combination of the set size of esco, with the intersection size with the others
set_data <- data.frame(
  'Hypomyces_Cladobotryum' = length(listInput$`Hypomyces_Cladobotryum`),
  Trichoderma = length(listInput$Trichoderma),
  Escovopsis = length(listInput$Escovopsis),
           Overlap = length(intersect(which(listInput$Escovopsis %in% listInput$`Hypomyces_Cladobotryum`), 
                                      which(listInput$Escovopsis %in% listInput$Trichoderma)))) %>%
  gather(set, n_genes)

set_data$set <- factor(set_data$set, levels =c('Trichoderma', 'Hypomyces_Cladobotryum', 'Escovopsis', 'Overlap'))

ggplot(set_data, aes(y = n_genes, x = set, fill = set)) + 
  geom_col() +
  theme_classic() + 
  scale_fill_manual(values = clade_colors) +
  labs(y = '', 
       x = '') +
  geom_text(aes(y = n_genes +80, label = n_genes)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, max(set_data$n_genes)+200)) +
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        legend.justification=c(1,1),
        ) +
  ggsave('plots/esco_v_overlap_barchart.pdf', width = 4.5, height = 4)













###############
## ADD in functional annotations to our shared orthogroups
#################


fais <- list.files('annotation/proteins/', pattern = 'fai', full.names = TRUE)
aa_lengths <- lapply(fais, function(x){
  read_tsv(x, col_names = c('gene', 'length', NA, NA, NA)) %>%
    select(-contains('X'))
}) %>% bind_rows()

groups = list(c('Higher', 'Leafcutter', 'Lower', 'AllOutgroup'),
              c('Lower','AllOutgroup'),
              c('Coral', 'Higher', 'Lower', 'AllOutgroup'),
              c('Coral', 'Higher', 'Leafcutter', 'AllOutgroup'),
              c('Coral', 'Leafcutter', 'Lower', 'AllOutgroup'),
              c('Coral', 'Lower', 'AllOutgroup'),
              c('Higher', 'Leafcutter', 'AllOutgroup'),
              c('Leafcutter'),
              c('Lower'),
              c('Higher', 'Leafcutter'),
              c('Higher', 'Lower', 'AllOutgroup'),
              c('Coral', 'Higher', 'Leafcutter', 'Lower'),
              c('Leafcutter', 'Lower', 'AllOutgroup'))

groups <- lapply(groups, function(x){
  paste(sort(x), collapse = ",")
}) %>% unlist()



input = listInput[c(ag, 'AllOutgroup')]
elements <- unique(unlist(input))


data <- unlist(lapply(input, function(x) {
  x <- as.vector(match(elements, x))
}))

data[is.na(data)] <- as.integer(0)
data[data != 0] <- as.integer(1)
data <- data.frame(matrix(data, ncol = length(input), byrow = F))
data <- data[which(rowSums(data) != 0), ]
names(data) <- paste0("pa_", names(input))
rownames(data) <- elements
data$Orthogroup <- elements




annot_add <- annot %>%
   left_join(., ortho_full, by = 'gene') %>%
   left_join(., aa_lengths, by = 'gene') %>%
   left_join(data, ., by = c('Orthogroup')) %>%
  group_by(Orthogroup)









counts <- c()
data <- data[, paste0("pa_", c(ag, 'AllOutgroup'))]
cols <- colnames(data[, paste0('pa_', c(ag, 'AllOutgroup'))])


for (x in rownames(data)) {
  id = sort(cols[which(data[x, ] == 1)])
  id = paste(id, collapse = ",") %>% gsub("pa_", "", .)
  if(!id %in% names(counts)){
    counts[id] = x
  }  else {
    counts[id] = paste(counts[id], x, sep = ",")
  }
}

counts <- counts[groups]

counts_df <- lapply(names(counts), function(x){
  genes = strsplit(counts[x], split = ',')[[1]]
  data.frame(group = rep(x, length(genes)), Orthogroup = genes)
}) %>% bind_rows()




orthogroup_meta <- annot_add %>% 
  group_by(Orthogroup, tool) %>%
  mutate(combined_annot = paste(unique(annot), collapse = ",")) %>%
  ungroup()  %>%
  group_by(Orthogroup) %>%
  mutate(average_length = mean(length)) %>%
  ungroup() %>%
  select(Orthogroup,
         tool, 
         combined_annot,
         contains('pa_'),
         average_length) %>%
  distinct() %>%
  left_join(., counts_df, by = 'Orthogroup')








#orthogroup_meta %>%  write_excel_csv(., path = 'tables/geneorthology_pa_upset_data.csv')

return_secreted <- function(.){
  ## signalp == 1
  ## TMHMM == 1
  ## WolfPsort ==1
  ## targetp == 1
  ## AND
  ## dbCAN, merops, LED, virulence ID
  ## AND <= 300 amino acids
  ## used in a groupby dplyr statment as a mutate
  flag = 0
  if(.$SignalP_EUK == 1 && !is.na(.$SignalP_EUK)){ flag = 1 }
  else if (.$TMHMM == 1 && !is.na(.$TMHMM)){ flag = 1 }
  else if (.$wolfpsort == 1 && !is.na(.$wolfpsort)) { flag = 1}
  else if (.$targetp == 1 && !is.na(.$targetp)) {flag = 1}
  
  if(flag == 1){
    if(.$LED == 1 && !is.na(.$LED)) { flag = 1}
    else if (.$MEROPS == 1 && !is.na(.$MEROPS)) {flag = 1}
    else if (.$cazy == 1&& !is.na(.$cazy)) {flag = 1}
    else if (.$virulence == 1 && !is.na(.$virulence)) {flag = 1}
    else{ flag = 0 }
  }
  if (flag == 1) {
    if(.$average_length >= 300) {flag = 1}
    else {flag = 0}}
  else { flag = 0}
  return(flag)
}




tools_interest <- c('CARD', 
                    'cazy', 
                    'led', 
                    'MEROPS', 
                    'SignalP_EUK', 
                    'virulence', 
                    'antismash',
                    'other')

ortho_subset <- orthogroup_meta %>% 
  filter(!is.na(group)) %>%
  group_by(group, tool) %>%
  mutate(tool_keep = case_when(
    tool %in% tools_interest ~ tool,
    TRUE ~ 'other'
  )) %>%
  ungroup() %>%
  mutate(annot_pa = 1)  %>%
  select(Orthogroup, annot_pa, tool_keep, group) %>%
  filter(tool_keep != 'other') %>%
  distinct() %>%
  spread(tool_keep, annot_pa) %>%
  dplyr::rename('BGC' = antismash,
         'Resistance' = CARD,
         'CAZyme' = cazy,
         'Lipase' = led,
         'Peptidase' = MEROPS,
         'Signaling' = SignalP_EUK, 
         'Virulence' = virulence)
  
ortho_subset[is.na(ortho_subset)] <- 0
  
ortho_subset %>%
  mutate(group = list(strsplit(group, split = ',')[[1]])) %>%
  unnest_longer(group) %>%
  write_csv('tables/orthogroup_group_functional_annotation.csv')


pca_data <- ortho_subset %>% 
  dplyr::select(-group) %>%
  distinct() %>%
  column_to_rownames(var = 'Orthogroup')
         



library(vegan)
pca_out    <- prcomp(pca_data)
eigs       <- pca_out$sdev^2
proportion <- (eigs/sum(eigs))*100
cumulative <- cumsum(eigs)/sum(eigs)
screeplot(pca_out)
n_k <- length(which(proportion >= 10))

if(n_k < 2){
  n_k = 2
}

## nmds
## bray curtis distance is more resilient to nulls
example_NMDS=metaMDS(pca_data, k = n_k) # The number of reduced dimensions ## components with >10% variance explained
stressplot(example_NMDS)

## test for significance
c_dist <- vegdist(pca_data)




#ano_test <- anosim(c_dist, grouping = treat)
#summary(ano_test)
#plot(ano_test)

#
groups_names <- unique(ortho_subset$group)
groups_n <- seq(1, length(groups_names))
names(groups_n) = groups_names

nmds <- data.frame(scores(example_NMDS)) %>%
  rownames_to_column( var = 'Orthogroup') %>%
  left_join(ortho_subset) %>%
  group_by(group) %>% 
  mutate(group_n = groups_n[group]) %>%
  ungroup()
  



p <- ggord::ggord(example_NMDS, ellipse = FALSE, size = 2, color = 'black') +
  theme_classic() +
  scale_shape_manual('Groups', values = seq(1, 13)) 

p$layers[[1]] <- NULL

## 76 unique combinations 
nmds %>% dplyr::select(-Orthogroup, -contains('NMDS'), -group, -group_n) %>% 
  distinct() %>%
  nrow()

double_dig <- nmds %>%
  filter(group_n >= 10)
single_dig <- nmds %>%
  filter(group_n < 10)

p + geom_text(data = single_dig, aes(x = NMDS1, y = NMDS2, group = group, label = group_n), size = 3, check_overlap = FALSE, nudge_x = 0.02) +
  geom_text(data = double_dig, aes(x = NMDS1, y = NMDS2, group = group, label = group_n), size = 3, check_overlap = FALSE, nudge_x = -0.025) +
  geom_point(alpha = 0.2)

  
nmds_map <- nmds %>% dplyr::select(-contains('NMDS'), -Orthogroup)

# nmds %>% 
#   dplyr::select(-Orthogroup, -NMDS3, -NMDS4, -group, group_n) %>%
#   mutate_at(vars(contains('NMDS')), signif, digits = 2) %>%
#   dplyr::group_by(antismash, CARD, cazy, led, MEROPS, SignalP_EUK, virulence) %>%
#   mutate(group_nest = paste(unique(group_n), collapse = ',')) %>%
#   ungroup() %>%
#   dplyr::select(-group_n) %>%
#   distinct() %>%
#   data.frame() %>% 
#   head()
  #%>%
#  left_join(., nmds_map, by = c('antismash', 'CARD', 'cazy', 'led', 'MEROPS', 'SignalP_EUK', 'virulence'))
    
nmds_regroup <- nmds %>%
  mutate(group = ifelse(group %in% c('Lower', 'Leafcutter'), yes = group, no = 'other'))

p + geom_jitter(data = nmds_regroup, alpha = 0.6, aes(x = NMDS1, y = NMDS2, colour = group), width = 0.3, height = 0.3, size = .2) +
  scale_colour_manual(values = c('Lower'= '#e6ab02', 'Leafcutter' = '#0069a3', 'other' = 'black')) +
  ggsave('plots/ordination_diverse_orthogroups_annotation.pdf', height = 5, width = 5)

p + geom_point(data = nmds_regroup, alpha = 0.8, aes(x = NMDS1, y = NMDS2, colour = group)) +
  scale_colour_manual(values = c('Lower'= '#e6ab02', 'Leafcutter' = '#0069a3', 'other' = 'black'))


nmds %>%
  dplyr::select(Orthogroup, group, group_n, contains('NMDS'), everything()) %>%
  write_csv(., 'tables/ordination_diverse_orthogroups_annotation_legend.csv')


# secreted <- ortho_subset %>%
#   rowwise() %>%
#   do(secreted = return_secreted(.)) %>% 
#   .$secreted %>% unlist()
#   
# 
# ortho_subset$secreted <- secreted
# 
# 
# ortho_subset <- ortho_subset %>%
 #  gather(tool, combined_annot, -OrthoGroup, -contains('pa_'), -group,)

# n_genes_data <- ortho_subset %>% select(OrthoGroup, group) %>%
#   group_by(group) %>%
#   summarize(n_genes = n_distinct(OrthoGroup)) %>%
#   mutate( y = 25)



# ortho_subset %>%
#   ggplot(.) +
#   geom_bar(stat = 'count', aes(x = group, fill = tool_keep)) +
#   scale_fill_brewer(palette = 'Dark2') +
#   theme_classic() + 
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.title = element_blank(),
#         legend.position = c(.9,.6))  +
#   scale_y_continuous(expand = c(0,0)) +
#   labs(y = '', x = '') +
#   ggsave('plots/orthologues_upset_subset_function_annot.pdf', width = 4, height = 5)

## grab LED, merops, cazy, secreted, signalp, virulence, tmhmm








## Analysis ends here! Nothing gets enriched


go_annot <- annot_add %>% filter(tool == 'GO') %>%  
  select(Orthogroup,
         annot,
         contains('pa_'), Agriculture) %>% 
  group_by(Orthogroup) %>%
  mutate(AllAgs = paste(sort(unique(Agriculture)), collapse = ',')) %>%
  ungroup() %>% select(-Agriculture) %>%
  distinct() %>%
  gather(agric, val, -Orthogroup, -annot, -AllAgs) %>%
  select(-agric) %>%
  group_by(Orthogroup) %>%
  mutate(go = paste(unique(strsplit(annot, ",")[[1]]), collapse = ",")) %>%
  ungroup() %>%
  select(-annot) %>%
  distinct()

go_annot_u      <- go_annot %>% select(go, Orthogroup) %>% distinct()
ortho2go        <- go_annot_u$go
names(ortho2go) <- go_annot_u$Orthogroup


ortho2go        <- lapply(ortho2go, function(x){strsplit(x, ',')[[1]]})
geneNames        <- names(ortho2go)



groups = list(c('Higher', 'Leafcutter', 'Lower', 'AllOutgroup'),
              c('Lower','AllOutgroup'),
              c('Coral', 'Higher', 'Lower', 'AllOutgroup'),
              c('Coral', 'Higher', 'Leafcutter', 'AllOutgroup'),
              c('Coral', 'Leafcutter', 'Lower', 'AllOutgroup'),
              c('Coral', 'Lower', 'AllOutgroup'),
              c('Higher', 'Leafcutter', 'AllOutgroup'),
              c('Leafcutter'),
              c('Lower'),
              c('Higher', 'Leafcutter'),
              c('Higher', 'Lower', 'AllOutgroup'),
              c('Coral', 'Higher', 'Leafcutter', 'Lower'),
              c('Leafcutter', 'Lower', 'AllOutgroup'))

groups <- lapply(groups, function(x){
  paste(sort(x), collapse = ",")
}) %>% unlist()


colnames(data) <- sub('pa_', '', colnames(data))
counts <- c()
data <- data[, c(ag, 'AllOutgroup')]
cols <- colnames(data[, c(ag, 'AllOutgroup')]) 


for (x in rownames(data)) {
  id = sort(cols[which(data[x, ] == 1)])
  id = paste(id, collapse = ",")
  if(!id %in% names(counts)){
    counts[id] = x
  }  else {
    counts[id] = paste(counts[id], x, sep = ",")
  }
}

xx <- as.list(GOTERM)
library(topGO)
go_enrich <- lapply(groups, function(x){
  ## define the genes that I'm interested in (those not in fungus farming ant escovopsis)
  myInterestingGenes <- strsplit(counts[x], split = ",")[[1]]
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
})




### KEGG enrichment test



# path <- keggList("pathway")
# head(path)
# 
# ko_to_path <- keggLink("ko", 'pathway')
# ko_path_df <- data.frame(KO = ko_to_path, path_id = names(ko_to_path))
# 
# 
# 
# keggframeData = distinct(data.frame(kegg = sub('ko:K', '', annot$KO), gene = annot$OrthoGroup, 
#                          stringsAsFactors = FALSE))
# head(keggframeData)
# keggFrame=KEGGFrame(keggframeData, organism = 'Escovopsis')
# 
# gsc <- GeneSetCollection(keggFrame, setType = KEGGCollection())
# 
# universe = unique(keggframeData$gene)
# 
# 
# genes = orthogroup_meta %>% filter(group == 'Leafcutter') %>% .$OrthoGroup
# 
# kparams <- GSEAKEGGHyperGParams(name="My Custom GSEA based annot Params",
#                                   geneSetCollection=gsc,
#                                   geneIds = genes,
#                                   universeGeneIds = universe,
#                                   pvalueCutoff = 0.05,
#                                   testDirection = "over")
# kOver <- hyperGTest(kparams)
# 
# results <- summary(kOver) %>%
#   data.frame() %>%
#   mutate(BH.correction = signif(p.adjust(Pvalue,method="BH"), 6))  %>%
#   mutate(KO = paste0('ko:K', KEGGID)) %>%
#   left_join(., ko_path_df, 'KO') %>%
#   mutate(path_def = path[path_id]) %>%
#   filter(!is.na(path_def),
#          BH.correction <= 0.05)
#    
# head(arrange(results, Pvalue), 50) %>% .$path_def

# venn.diagram(
#   x = list('Escovopsis' = esco_set, 'Hypomyces/Cladobotryum' = ch_set, 'Trichoderma' = t_set),
#   height = 3000,
#   width = 3000,
#    category.names = c(paste0("Escovopsis\n(", n_distinct(esco_set), " Orthologues across ", nrow(filter(meta, clade_groups == 'Escovopsis')), " genomes)"),
#                       paste0("Cladobotryum_Hypomyces\n(", n_distinct(ch_set), " Orthogroups across ", nrow(filter(meta, clade_groups == 'Hypomyces/Cladobotryum')), " genomes)"), 
#                       paste0("Trichoderma\n(", n_distinct(t_set), " Orthogroups across ", nrow(filter(meta, clade_groups == 'Trichoderma')), " genomes)")),
#   filename = 'plots/shared_genes_venn_diagramm.png',
#   resolution = 300,
#   # Numbers
#   fontface = "bold",
#   fontfamily = "sans",
#   main.fontface = 'bold',
#   main.fontfamily = 'sans',
#   # Set names
#   cat.fontface = "bold",
#   cat.default.pos = "outer",
#   cat.pos = c(-27, 27, 135),
#   cat.dist = c(0.055, 0.055, 0.085),
#   cat.fontfamily = "sans",
#   rotation = 1
# )



# ortho_annot <- ortho_full %>%
#   mutate(cazy = ifelse(gene %in% cazy,
#                        yes = 'cazy',
#                        no = NA)) %>%
#   mutate(card = ifelse(gene %in% card,
#                        yes = 'card',
#                        no = NA)) %>%
#   select(-gene, -genome) %>%
#   filter(!is.na(cazy) | !is.na(card)) %>%
#   group_by(OrthoGroup) %>%
#   mutate(annotation = paste(na.omit(c(unique(cazy), unique(card))), collapse = "_")) %>%
#   ungroup()  %>%
#   select(-cazy, -card) %>%
#   distinct() %>%
#   data.frame(., stringsAsFactors = FALSE)
#   
# 
# 
# 
# 
# 
# rownames(ortho_full) <- ortho_full$OrthoGroup
# ortho_full$OrthoGroup <- NULL
# 
# 
# ortho <- read_csv('annotation/fastortho/orthologues_presence_absence.csv')
# ortho_mat <- data.frame(ortho, stringsAsFactors = FALSE)
# rownames(ortho_mat) <- ortho$OrthoGroup
# ortho_mat$OrthoGroup <- NULL
# ortho_mat <- as.matrix(ortho_mat)
# 
# 
# ortho_mat[ortho_mat > 1] <- 1
# 
# 
# breaks = seq(0, 1)
# 
# #colors_f <- colorRampPalette(colors = c('white', 'black'))
# colors <- c('white', 'black')
# 
# 
# annot_colors <- list(annotation = c('card' =  "#1B9E77", 'cazy' = "#D95F02", "cazy_card" = "#7570B3"))
# 
# pheatmap(ortho_mat, show_rownames = FALSE, color = colors, cellwidth = 10, annotation_row = ortho_full, annotation_colors = annot_colors, 
#          annotation_names_row = FALSE, filename = 'plots/orthologous_gene_matrix.pdf')
# 
# 
