library(KEGGREST)
library(GOstats)
library(GSEABase)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)
library(UpSetR)


source('scripts/color_palettes.R')

annot <- read_tsv('annotation/all_annotations.txt', col_names = c('acc', 'gene', 'tool', 'annot'))

>>>>>>> 602f3b502c4c24b642d4823f1aa7308a586139bd

meta <- read_csv('tables/metadata.csv')

all_gene_id_files <- list.files('annotation/fastortho/proteins/', pattern = 'gene_ids', full.names = TRUE)

all_gene_ids <- lapply(all_gene_id_files, function(x){
  ids <- scan(x, what = 'character')
  data.frame(gene = ids, genome = rep(basename(x), length(ids)), stringsAsFactors = FALSE) %>%
    mutate(acc = case_when(
      grepl('GCA_004303015', genome) ~ 'GCA_004303015.1',
      grepl('GCA', genome) ~ paste(strsplit(genome, split = "_")[[1]][c(1,2)], collapse = "_"),
      grepl("SPDT00000000.1_genomic", genome) ~ "SPDT00000000",
      grepl('LGSR', genome) ~ "LGSR00000000",
      grepl('.fasta_gene_ids', genome) ~ sub(".all.maker.proteins.fasta_gene_ids", '', genome),
      TRUE ~ genome
    )) %>%
    select(-genome) 
}) %>% bind_rows()

ortho_full <- read_csv('annotation/fastortho/orthologues_full.csv', col_names = c('OrthoGroup', 'genome', 'gene')) %>%
  group_by(genome) %>%
  mutate(acc = case_when(
    grepl('GCA_004303015', genome) ~ 'GCA_004303015.1',
    grepl('GCA', genome) ~ paste(strsplit(genome, split = "_")[[1]][c(1,2)], collapse = "_"),
    grepl("SPDT00000000.1_genomic", genome) ~ "SPDT00000000",
    grepl('all.maker', genome) ~ sub(".all.maker.proteins", '', genome),
    grepl('LGSR', genome) ~ "LGSR00000000",
    TRUE ~ genome
  )) %>%
  ungroup() %>%
  select(-genome) %>%
  full_join(all_gene_ids, ., by = c('gene', 'acc')) %>%
  left_join(., meta, by = 'acc') %>%
  mutate(orthology_meta = ifelse(is.na(OrthoGroup),
                                 yes = gene,
                                 no = OrthoGroup))


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
    group_by(orthology_meta) %>%
    mutate(n_genomes = length(unique(genus_species))) %>%
    ungroup() %>%
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


# esco_set <- venn %>% filter(clade_groups == 'Escovopsis') %>% .$orthology_meta %>% unique()
# ch_set   <- venn %>% filter(clade_groups == 'Hypomyces/Cladobotryum') %>% .$orthology_meta %>% unique()
# t_set    <- venn %>% filter(clade_groups == 'Trichoderma') %>% .$orthology_meta %>% unique()

listInput <- lapply(c(groups, ag), 
                    function(x){
                      venn %>% filter(clade_groups == x) %>% .$orthology_meta %>% unique()
                    })

names(listInput) <- c(groups, ag)

listInput[['AllOutgroup']] <- c(listInput$Trichoderma, listInput$`Hypomyces/Cladobotryum`)

# listInput <- list('Escovopsis' = esco_set,
#                   'Hypomyces/Cladobotryum' = ch_set,
#                   'Trichoderma' = t_set)
# set_name(listInput)

pdf("plots/orthologues_upset_summary.pdf", width = 8, height = 8)
upset(fromList(listInput), order.by = "freq", 
      sets = c(groups),
       point.size = 3.5, 
       line.size = 1.5, 
       text.scale = c(1.3, 1.3, 1, 1, 2, 0.75),
       mb.ratio = c(0.8, 0.2))
dev.off()

pdf("plots/orthologues_upset_subset.pdf", width = 8, height = 8)

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
  'Hypomyces.Cladobotryum' = length(listInput$`Hypomyces/Cladobotryum`),
  Trichoderma = length(listInput$Trichoderma),
  Escovopsis = length(listInput$Escovopsis),
           Overlap = length(intersect(which(listInput$Escovopsis %in% listInput$`Hypomyces/Cladobotryum`), 
                                      which(listInput$Escovopsis %in% listInput$Trichoderma)))) %>%
  gather(set, n_genes)

set_data$set <- factor(set_data$set, levels =c('Trichoderma', 'Hypomyces.Cladobotryum', 'Escovopsis', 'Overlap'))

ggplot(set_data, aes(y = n_genes, x = set)) + geom_col(fill = 'black') +
  theme_minimal() + 
  labs(y = 'Number of Orthologous Genes per Clade', 
       x = '') +
  ggsave('plots/esco_v_overlap_barchart.pdf', width = 5, height = 7)


venn.diagram(
  x = list('Escovopsis' = esco_set, 'Hypomyces/Cladobotryum' = ch_set, 'Trichoderma' = t_set),
  height = 3000,
  width = 3000,
   category.names = c(paste0("Escovopsis\n(", n_distinct(esco_set), " Orthologues across ", nrow(filter(meta, clade_groups == 'Escovopsis')), " genomes)"),
                      paste0("Cladobotryum_Hypomyces\n(", n_distinct(ch_set), " Orthogroups across ", nrow(filter(meta, clade_groups == 'Hypomyces/Cladobotryum')), " genomes)"), 
                      paste0("Trichoderma\n(", n_distinct(t_set), " Orthogroups across ", nrow(filter(meta, clade_groups == 'Trichoderma')), " genomes)")),
  filename = 'plots/shared_genes_venn_diagramm.png',
  resolution = 300,
  # Numbers
  fontface = "bold",
  fontfamily = "sans",
  main.fontface = 'bold',
  main.fontfamily = 'sans',
  # Set names
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)
=======
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

listInput[['AllOutgroup']] <- c(listInput$Trichoderma, listInput$`Hypomyces_Cladobotryum`)


pdf("plots/orthologues_upset_summary.pdf", width = 5, height = 4.5)
upset(fromList(listInput), order.by = "freq", 
      sets = c(groups),
       point.size = 3.5, 
       line.size = 1.5, 
       text.scale = c(1.3, 1.3, 1, 1, 2, 0.75),
      mb.ratio = c(0.8, 0.2))
dev.off()



pdf("plots/orthologues_upset_subset.pdf", width = 5, height = 6)
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
  ggsave('plots/esco_v_overlap_barchart.pdf', width = 3.5, height = 2)







## output a table with orthologous gene metadata


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
data$OrthoGroup <- elements


fais <- list.files('annotation/fastortho/proteins/', pattern = 'fai', full.names = TRUE)

aa_lengths <- lapply(fais, function(x){
  read_tsv(x, col_names = c('gene', 'length', NA, NA, NA)) %>%
    select(-contains('X'))
}) %>% bind_rows()





annot <- annot %>%
  left_join(., ortho_full, by = 'gene') %>% 
  left_join(., aa_lengths, by = 'gene') %>%
  left_join(data, ., by = c('OrthoGroup'))


# annot <- read_tsv('annotation/all_annotations.txt', col_names = c('genome', 'gene', 'tool', 'annot')) %>%
#   left_join(., ortho_full, by = 'gene') %>% 
#   left_join(., aa_lengths, by = 'gene') %>%
#   left_join(data, ., by = c('OrthoGroup'))


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
  data.frame(group = rep(x, length(genes)), OrthoGroup = genes)
}) %>% bind_rows()




orthogroup_meta <- annot %>% 
  group_by(OrthoGroup, tool) %>%
  mutate(combined_annot = paste(unique(annot), collapse = ",")) %>%
  ungroup()  %>%
  group_by(OrthoGroup) %>%
  mutate(average_length = mean(length)) %>%
  ungroup() %>%
  select(OrthoGroup,
         tool, 
         combined_annot,
         contains('pa_'),
         average_length) %>%
  distinct() %>%
  left_join(., counts_df, by = 'OrthoGroup')

orthogroup_meta %>%  write_excel_csv(., path = 'tables/geneorthology_pa_upset_data.csv')

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
                    'LED', 
                    'MEROPS', 
                    'SignalP_EUK', 
                    'virulence', 
                    'other')
ortho_subset <- orthogroup_meta %>% 
  filter(!is.na(group)) %>%
  group_by(group, tool) %>%
  mutate(tool_keep = case_when(
    tool %in% tools_interest ~ tool,
    TRUE ~ 'other'
  )) %>%
  ungroup() %>%
  mutate(annot_pa = 1)
  


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

n_genes_data <- ortho_subset %>% select(OrthoGroup, group) %>%
  group_by(group) %>%
  summarize(n_genes = n_distinct(OrthoGroup)) %>%
  mutate( y = 25)



ortho_subset$group <- factor(ortho_subset$group, levels = groups)



ortho_subset %>%
  select(OrthoGroup, annot_pa, tool_keep, group) %>%
  filter(tool_keep != 'other') %>%
  distinct() %>%
  ggplot(.) +
  geom_bar(stat = 'count', aes(x = group, fill = tool_keep)) +
  scale_fill_brewer(palette = 'Dark2') +
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(.9,.6))  +
  scale_y_continuous(expand = c(0,0)) +
  labs(y = '', x = '') +
  ggsave('plots/orthologues_upset_subset_function_annot.pdf', width = 4, height = 5)

## grab LED, merops, cazy, secreted, signalp, virulence, tmhmm








## Analysis ends here! 
### Just keeping this here because only one group was enriched:: 

## Leafcutters for:

# go_id                Term Ontol        pval BH.correction
# GO:0007155 GO:0007155       cell adhesion    BP 5.35002e-06    0.00558007
# GO:0022610 GO:0022610 biological adhesion    BP 5.35002e-06    0.00558007


go_annot <- annot %>% filter(tool == 'GO') %>%  
  select(OrthoGroup,
         annot,
         contains('pa_'), 
         AllAgs) %>% 
  gather(agric, val, -OrthoGroup, -annot, -AllAgs) %>%
  select(-agric) %>%
  group_by(OrthoGroup) %>%
  mutate(go = paste(unique(strsplit(annot, ",")[[1]]), collapse = ",")) %>%
  ungroup() %>%
  select(-annot) %>%
  distinct()

go_annot_u      <- go_annot %>% select(go, OrthoGroup) %>% distinct()
ortho2go        <- go_annot_u$go
names(ortho2go) <- go_annot_u$OrthoGroup


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



path <- keggList("pathway")
head(path)

ko_to_path <- keggLink("ko", 'pathway')
ko_path_df <- data.frame(KO = ko_to_path, path_id = names(ko_to_path))

kegg_f <- list.files('annotation/kofamscan', full.names = TRUE, pattern = 'fasta_fix')

kegg <- lapply(kegg_f, function(x){
  read_delim(x, col_names = c('gene', 'KO', NA, NA, 'eval', 'def'), skip = 2, delim = " ", 
             trim_ws = TRUE,
             col_types = 'ccccnc') %>%
    select(-X3, -X4) %>%
    mutate(genome = basename(x)) %>% 
    filter(KO != "-")
}) %>% bind_rows() %>%
  mutate(genome = case_when(
    'LGSR' %in% genome ~ 'LGSR00000000',
    TRUE ~ sub('.all.maker.proteins.fasta_fix', '', genome)
  )) %>%
  mutate(KO = paste0('ko:', KO)) %>%
  left_join(., ko_path_df, by = 'KO') %>%
  mutate(path_def = path[path_id]) %>%
  filter(!is.na(path_def))

annot <- kegg %>%
  left_join(., ortho_full, by = 'gene') %>% 
  left_join(., aa_lengths, by = 'gene') %>%
  left_join(data, ., by = c('OrthoGroup'))


orthogroup_meta <- annot %>% 
  group_by(OrthoGroup) %>%
  mutate(average_length = mean(length)) %>%
  ungroup() %>%
  dplyr::select(OrthoGroup,
         contains('pa_'),
         average_length) %>%
  distinct() %>%
  left_join(., counts_df, by = 'OrthoGroup')

keggframeData = distinct(data.frame(kegg = sub('ko:K', '', annot$KO), gene = annot$OrthoGroup, 
                         stringsAsFactors = FALSE))
head(keggframeData)
keggFrame=KEGGFrame(keggframeData, organism = 'Escovopsis')

gsc <- GeneSetCollection(keggFrame, setType = KEGGCollection())

universe = unique(keggframeData$gene)


genes = orthogroup_meta %>% filter(group == 'Leafcutter') %>% .$OrthoGroup

kparams <- GSEAKEGGHyperGParams(name="My Custom GSEA based annot Params",
                                  geneSetCollection=gsc,
                                  geneIds = genes,
                                  universeGeneIds = universe,
                                  pvalueCutoff = 0.05,
                                  testDirection = "over")
kOver <- hyperGTest(kparams)

results <- summary(kOver) %>%
  data.frame() %>%
  mutate(BH.correction = signif(p.adjust(Pvalue,method="BH"), 6))  %>%
  mutate(KO = paste0('ko:K', KEGGID)) %>%
  left_join(., ko_path_df, 'KO') %>%
  mutate(path_def = path[path_id]) %>%
  filter(!is.na(path_def),
         BH.correction <= 0.05)
   
head(arrange(results, Pvalue), 50) %>% .$path_def

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
>>>>>>> 602f3b502c4c24b642d4823f1aa7308a586139bd



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
