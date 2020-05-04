library(tidyverse)
library(cowplot)

source('scripts/color_palettes.R')
metadata <- read_csv('tables/metadata.csv') %>% data.frame()



proteins <- dplyr::filter(metadata, !is.na(n_proteins)) %>%
  ggplot(aes(x = clade_groups, y = n_proteins, colour = clade_groups)) + 
  geom_boxplot() + 
  theme_classic() +
  labs(x = '', 
       y = 'Number of Genes')  +
  scale_y_continuous(limits = c(0, max(metadata$n_proteins))) +
  theme(
        legend.position = 'none',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = clade_colors) + 
  ggsave('plots/protein_coding_gene_N_boxplot.pdf', width = 3, height = 3)

as_len <- metadata %>% filter(!is.na(assembly_length)) %>%
  ggplot(aes(x = clade_groups, y = assembly_length, colour = clade_groups)) + 
  geom_boxplot() + 
  theme_classic() +
  labs(x = '', 
       y = 'Assembly Length')  +
  scale_y_continuous(limits = c(0, max(metadata$assembly_length))) +
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  scale_colour_manual(values = clade_colors) +
  ggsave('plots/assembly_length_boxplot.pdf', width = 3, height = 3)


gs_shape_data <- filter(metadata, 
                        GenomeSizeEstimationStrategy %in% c('FI-PC', 'pacbio'))
genome_size <- 
  metadata %>% filter(!is.na(GenomeSize)) %>%
  ggplot(aes(x = clade_groups, 
             y = GenomeSize, 
             colour = clade_groups)) + 
  geom_boxplot() + 
  theme_classic() +
  scale_colour_manual(values = clade_colors) +
  geom_jitter(data = gs_shape_data, aes(shape = GenomeSizeEstimationStrategy),
              colour = 'black', width = 0.25) +
  labs(x = '', 
       y = 'Genome Size') +
  scale_y_continuous(limits = c(0, max(metadata$GenomeSize))) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  ggsave('plots/genome_size_boxplot.pdf', width = 10, height = 3)
genome_size

plot_grid(as_len, genome_size, labels = c('A', 'B'), label_size = 12) + 
  ggsave2(filename = 'plots/supplementary_fig1.pdf')



### stats files were made with agat with the following command
### agat_sp_statistics.pl -o GCA_011799845.1.stats -gff /home/gotting/gene_annotation/GCA_011799845.1.genemodels.gff
### cat GCA_011799845.1.stats | sed -e"s/ /_/g" | sed -e "s/___*/ /" >GCA_011799845.1.stats.fix

gff_f <- list.files('gff_stats', full.names = TRUE, pattern = 'fix')
gff_stats <- lapply(gff_f, function(x){
  read_delim(x, delim = " ", col_names = c('stat','value')) %>% 
    mutate(filename = sub("_genomic", '', sub(".stats.fix", "", basename(x))))
}) %>% bind_rows()
gff_stats$acc <- sapply(gff_stats$filename, function(x){
  metadata[grep(x, metadata$genome_id), 'acc']
})


gff_stats <- gff_stats %>% left_join(., metadata, by = 'acc')

stats_interest <- c("mean_cds_length" , 'Total_cds_length', 'Number_of_single_exon_gene')

sub <- gff_stats %>% filter(stat %in% stats_interest)
sub$stat <- factor(sub$stat, levels = stats_interest)

sub %>% filter(stat %in% stats_interest) %>%
  ggplot(aes(x = clade_groups, y = value, colour = clade_groups)) + 
  geom_boxplot() + 
  theme_classic() +
  facet_wrap(~stat, scales = 'free') +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  scale_colour_manual(values = clade_colors) +
  ggsave('plots/lengths_boxplot.pdf')

stats_interest <- c('Number_of_intron_in_cds',
                    "mean_introns_in_cdss_per_mrna",
                    'mean_intron_in_cds_length',
                    "Total_intron_length_per_cds"
                    )



sub <- gff_stats %>% filter(stat %in% stats_interest)
sub$stat <- factor(sub$stat, levels = stats_interest)

ggplot(sub, aes(x = clade_groups, y = value, colour = clade_groups)) + 
  geom_boxplot() + 
  theme_classic() +
  facet_wrap(~stat, scales = 'free') +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  scale_colour_manual(values = clade_colors) +
  ggsave('plots/introns_boxplot.pdf')

##source('scripts/make_map.R')
##subplot(map, proteins, proteins, proteins, nrows = 2) 



