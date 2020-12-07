library(tidyverse)
library(ggtree)
library(treeio)

load('rdata/orthologues_long_annot.rdata')

head(orthog_long_annot)

ani <- read_tsv('tables/ANIvis.tsv') %>% 
  mutate_at(c('Query', 'Reference'), basename) %>%
  mutate_at(c('Query', 'Reference'), function(x){sub(".cds_from_genomic.fa", "", x)}) %>%
  mutate_at(c('Query', 'Reference'), function(x){sub("\\..*$", "", x)})



tree <- read.tree("astral.nwk")
tips <- tree$tip.label
names(tips) <- sub("\\..*$", '', tips)

genomes <- orthog_long_annot$acc %>% unique()

tips_remove = tips[!names(tips) %in% genomes]

genomes_keep <- genomes[genomes %in% names(tips)]

tree_sub <- drop.tip(tree, tips_remove)

ortho <- orthog_long_annot %>% 

  group_by(Orthogroup, acc) %>%
  mutate(genome_copies = length(unique(gene))) %>%
  select(-genome, -gene) %>%
  filter(tool %in% c('CARD', 'MEROPS', 'cazy',
                     'led', 'antismash', 'busco'),
         !grepl("mix", annot)) %>%
  distinct() %>%
  group_by(Orthogroup) %>% 
  mutate(num_annots = length(unique(annot))) %>%
  ungroup()


load('rdata/gene_network.rdata')  
  

w

  
ani_sub <- ani %>%
  select(-QueryAligned) %>%
  filter(Reference %in% genomes_keep & Query %in% genomes_keep) %>% spread(Query, ANI) %>%
  column_to_rownames(var = 'Reference') 
ani_sub <- ani_sub[,tree_sub$tip.label[tree_sub$tip.label %in% colnames(ani_sub)]]


p <- ggtree(tree_sub) +
  theme_tree() +
  geom_tiplab() 


gene_data <- gene_network %>%
  ungroup() %>%
  mutate(genome = sub(".all.maker.proteins", '', genome)) %>%
  mutate(genome = sub("\\..*$", "", genome)) %>%
  filter(genome %in% genomes_keep) %>%
  mutate(num_genes = ifelse(num_genes > 3, yes = 3, no = num_genes)) %>%
  group_by(Orthogroup) %>%
  mutate(num_genomes = length(unique(genome))) %>%
  ungroup()


ggplot(head(gene_data, 1000), aes(y = genome, x = Orthogroup, fill = num_genes)) +
  geom_tile() +
  theme_classic() +
  coord_polar() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_blank())  #+
  #ggsave("test_heat.png")

gheatmap(p, 
         gene_data,
         offset=5, width=0.5, font.size=3, 
         colnames_angle=-45, hjust=0, 
         colnames = FALSE)









library(circlize)
## sample code



set.seed(123)
mat1 = rbind(cbind(matrix(rnorm(50*5, mean = 1), nr = 50), 
                   matrix(rnorm(50*5, mean = -1), nr = 50)),
             cbind(matrix(rnorm(50*5, mean = -1), nr = 50), 
                   matrix(rnorm(50*5, mean = 1), nr = 50))
)
rownames(mat1) = paste0("R", 1:100)
colnames(mat1) = paste0("C", 1:10)
mat1 = mat1[sample(100, 100), ] # randomly permute rows
split = sample(letters[1:5], 100, replace = TRUE)
split = factor(split, levels = letters[1:5])

ss


col_fun1 = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
circos.heatmap(mat1, split = split, col = col_fun1, show.sector.labels = TRUE)

circos.clear()





factors = letters[1:4]
circos.par("canvas.xlim" = c(-1, 1.5), "canvas.ylim" = c(-1, 1.5), start.degree = -45, points.overflow.warning = FALSE)
circos.initialize(factors = factors, xlim = c(0, 1))
circos.trackPlotRegion(ylim = c(0, 1), bg.col = NA, bg.border = NA)
circos.updatePlotRegion(sector.index = "a")
circos.heatmap(mat1, col = col_fun1)



circos.updatePlotRegion(sector.index = "b")
circos.text(0.5, 0.5, "first one", niceFacing = TRUE)

circos.clear()

par(new = TRUE)
circos.par("canvas.xlim" = c(-1.5, 1), "canvas.ylim" = c(-1.5, 1), start.degree = -45)
circos.initialize(factors = factors, xlim = c(0, 1))
circos.trackPlotRegion(ylim = c(0, 1), bg.col = NA, bg.border = NA)
circos.updatePlotRegion(sector.index = "d")


circos.text(0.5, 0.5, "second one", niceFacing = TRUE)
circos.updatePlotRegion(sector.index = "c")
circos.text(0.5, 0.5, "second one", niceFacing = TRUE)

circos.clear()

par(op)
