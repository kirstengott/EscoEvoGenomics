#!/usr/bin/env Rscript




library(tidyverse)

source(here::here('BGC_similarity_viz/make_BGC_plot.R'))


bgc_outdir <- 'BGC_plots'

if(!dir.exists(bgc_outdir)){
  dir.create(bgc_outdir)
}

## make basic plots of the BGCs for each component of the network
gff_files <- list.files('gff3', full.names = TRUE)

lapply(gff_files, function(x){
  plot_out <- paste0(bgc_outdir, '/', sub(".gff3", "", basename(x)), ".pdf")
  make_BGC_plot(gff3_file = x) +
    ggsave(plot_out, width = 40, height = 40, units = 'cm')
})