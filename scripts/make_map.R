
library(plotly)
library(tidyverse)
library(spData)


metadata <- readr::read_csv('tables/metadata.csv') %>%
  filter(genus == 'Escovopsis') %>%
  mutate(coordinates = sub(" ", "", coordinates)) %>%
  separate(coordinates, into = c('lat', 'long'), sep = ",") %>%
  mutate_at(c('lat', 'long'), as.numeric) %>%
  group_by(Country, Area) %>%
  mutate(size = length(unique(acc))) %>% #multiply by 5 for plotting
  ungroup()


south_america        = world[world$continent == "South America", ]
central_america_zone = world[world$name_long %in% c('Panama', 'Costa Rica'),]

total = rbind(south_america, central_america_zone) %>%
  rename('Country' = name_long) %>%
  left_join(., metadata, by = 'Country')

marker_sub <- total %>%
  filter(genus == 'Escovopsis')

g <- list(
  showland = TRUE,
  landcolor = toRGB("gray85"),
  scope = 'world',
  showcountries = TRUE,
  showrivers = TRUE,
  showlakes = TRUE,
  riverwidth = 3,
  rivercolor = '#0079B4',
  countrywidth = 2,
  countrycolor = '#636363',
  coastlinewidth = 2
)



map <- plot_geo(marker_sub,
        stroke = I("black"),
        type = "scattergeo",
        showlegend = FALSE,
        x = ~long, 
        y = ~lat,
        mode = 'markers',
        opacity = 0.5,
        marker = list(color = '#DB162F',
                      size = ~size*3,
                      showlegend = TRUE,
                      
                      line = list(
                        color = '#950f20',
                        width = 2
                      ))
) %>%
  layout(geo = g)

map

orca(map, "plots/map_plot.pdf")


















