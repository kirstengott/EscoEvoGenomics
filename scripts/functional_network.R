library(networkD3)
load('rdata/functional_network.rdata')
load('rdata/gene_network.rdata')

?forceNetwork


## nodes: annots
## links: orthogroups
## size: number of genomes
## group: detected in all, one, only escovopsis

links_map_n <- functional_network$Orthogroup %>% unique()

links_df <- lapply(links_map_n, function(x){
  annots <- functional_network %>% filter(Orthogroup == x) %>% .$annot %>% unique
  if(length(annots) < 2){
    data.frame(source = annots, target = annots)
  }
  else{
    links_ls <- combn(annots, m = 2, simplify = FALSE)
    lapply(links_ls, function(x){df = t(cbind(x))
      colnames(df) <- c('source', 'target')
      data.frame(df)
  }) %>% bind_rows()
  }
})

links_df <- bind_rows(links_df)




nodes <- functional_network %>% 
  select(annot, num_genomes, num_ortho) %>%
  distinct()


test_nodes <- nodes %>% head(., n = 10)

test_links <- links_df %>% filter(source %in% test_nodes$annot | target %in% test_nodes$annot) %>%
  mutate(value = 1)

forceNetwork(Links = test_links, Nodes = test_nodes,
             Source = "source", Target = "target",
             NodeID = "annot",
             Value = 'value',
             Group = "num_genomes", opacity = 0.8)

simpleNetwork(links_df)

data(MisLinks)
data(MisNodes)

forceNetwork(Links = MisLinks, Nodes = MisNodes,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",
             Group = "group", opacity = 0.8)

