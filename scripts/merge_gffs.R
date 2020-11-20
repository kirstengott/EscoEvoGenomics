#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
if (length(args) == 0L || any(c('-h', '--help') %in% args)) {
  message('usage: path/to/merge_gffs.R gff1 gff2 outfilename
    -h, --help        to print help messages')
  q('no')
}

library(rtracklayer)
library(GenomicRanges)
library(tidyverse)


gff1 <- rtracklayer::import(args[1])
gff2 <- rtracklayer::import(args[2])


gff_all <- lapply(unique(gff1$type), function(x){
  gff1_t = gff1[which(gff1$type == x), ]
  gff2_t = gff2[which(gff2$type == x), ]
  mergeByOverlaps(gff1_t, gff2_t)
})

gff <- do.call(rbind, gff_all)

data.frame(gff[, c('gff1_t')] )

save(gff, file = args[3])