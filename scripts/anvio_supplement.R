library('tidyverse')
source('scripts/color_palettes.R')

genes <- read_tsv('anvio/data.txt') %>% 
  gather(genome, presence, -Orthogroup) %>%
  filter(presence == 1)

data <- read_tsv('anvio/additional-items-data.txt') %>%
  gather(annotation, value, -Orthogroup) %>%
  filter(value == TRUE) %>%
  left_join(genes, data, by = "Orthogroup")






data %>%
  filter(!is.na(genome)) %>%
  ggplot(., aes(x = genome, y = presence)) +
  geom_col() +
  facet_grid(rows = vars(annotation), scales = 'free_y') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
