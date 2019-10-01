library(tidyverse)
library(pheatmap)




munge_genome_ids <- function(x){
  if(grepl('GCA', x)){
    paste(strsplit(x, split = "_")[[1]][c(1,2)], collapse = "_")
  } else {
    paste(strsplit(x, split = "_")[[1]][c(1)], collapse = "_")
  }
}





all_bigscape_f <- list.files(".", pattern = 'absence_presence')



all_input <- lapply(all_bigscape_f, function(x){
  type <- strsplit(x, "_")[[1]][1]
  read_tsv(x) %>%
  mutate(Type = type) %>%
  gather(FamilyID, Presence, -ACC, -Type)  %>%
  group_by(ACC) %>%
  mutate(Genome =  munge_genome_ids(ACC)) %>%
  ungroup()
})

# nrps <- read_tsv('./NRPS_absence_presence.tsv') %>%
#   mutate(Type = "NRPS") %>%
#   gather(FamilyID, Presence, -ACC, -Type)  %>%
#   group_by(ACC) %>%
#   mutate(Genome =  munge_genome_ids(ACC)) %>%
#   ungroup()
# 
# 
# pks <- read_tsv('PKSI_absence_presence.tsv') %>%
#   mutate(Type = 'PKS1')  %>%
#   gather(FamilyID, Presence, -ACC, -Type)  %>%
#   group_by(ACC) %>%
#   mutate(Genome =  munge_genome_ids(ACC)) %>%
#   ungroup()
# 
# terp <- read_tsv('terpene_absence_presence.tsv') %>%
#   mutate(Type = 'Terpene') %>%
#   group_by(ACC) %>%
#   mutate(Genome =  munge_genome_ids(ACC)) %>%
#   ungroup() %>%
#   gather(FamilyID, Presence, -ACC, -Type, -Genome)
# 





all <- bind_rows(all_input) %>%
  mutate(uniq_id = paste(Type, FamilyID, sep = ".")) %>%
  mutate(Presence = as.numeric(Presence)) %>%
  group_by(Type, FamilyID) %>%
  mutate(n_genomes = sum(Presence)) %>%
  ungroup() %>%
  filter(n_genomes > 2) %>%
  mutate(Presence = ifelse(Presence == 2, yes = 1, no = Presence))#%>% ## remove singletons
#  mutate(Genome = sub(pattern = 'Attine', replacement = 'Attini', Genome))

all %>% select(Type, FamilyID, Presence, Genome) %>% filter(Presence == 1) %>%
  group_by(Genome, Type) %>% 
  summarize(total = length(Type)) %>% write_csv(., path = "../../annotation/bigscape.csv")


heat <- all %>% select(Genome, uniq_id, Presence) %>%
  spread(uniq_id, Presence) %>%
  data.frame()

#heat[is.na(heat)] <- 0
rownames(heat) <- heat$Genome
heat$Genome <- NULL



metadata <- read_csv('metadata.csv') 

annotation_col = data.frame(
  "BGC_Type" = sub("\\..*$", "", colnames(heat_mat)),
  row.names = colnames(heat_mat)
)

annotation_row = data.frame(
 "Agriculture" = metadata$Agriculture,
row.names = rownames(heat_mat)
)

heat_mat <- as.matrix(heat[metadata$acc,])
rownames(heat_mat) <- metadata$genus_species

ann_colors = list(
'BGC_Type' = c(NRPS = '#1B9E77', PKS1 = '#D95F02', Terpene = '#7570B3'),
'Agriculture' = c(Leafcutter = 'green4', Lower = 'yellow', Coral = 'magenta3', Higher = 'blue3', 'NA' = 'white')
)

pheatmap(heat_mat, 
         annotation_col = annotation_col,
#         cluster_rows = FALSE,
        annotation_row = annotation_row,
         border_color = "grey70",
         cellwidth = 7,
         cellheight = 10,
         show_colnames = FALSE,
         annotation_colors = ann_colors,
         color = c('grey87', 'black'), legend = FALSE,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         filename = 'bgc_heat.pdf'
         #labels_row = sub("_.*$", "", rownames(heat_mat))
)




