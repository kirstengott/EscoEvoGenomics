library(tidyverse)

annot <- read_tsv('annotation/all_annotations.txt', col_names = c('acc', 'gene', 'tool', 'annot'))



tools_interest <- c('CARD', 
                    'cazy', 
                    'led', 
                    'MEROPS', 
                    'SignalP_EUK', 
                    'virulence', 
                    'antismash')

ref = 2

if(ref == 2){
  ref_col = 'g2'
  q_col = 'g1'
} else {
  ref_col = 'g1'
  q_col = 'g2'
}

mcscan_blocks <- read_tsv('mcscan/workingdir/ICBG742_NIGB00000000.i1.blocks', col_names = c('g1', 'g2')) %>%
  rename("gene" = all_of(ref_col)) %>%
  rename("gene_q" = all_of(q_col))

ref_bed <- read_tsv('mcscan/workingdir/NIGB00000000.bed', col_names = c('chr_r', 'start_r', 'end_r', 'gene'))
q_bed   <- read_tsv('mcscan/workingdir/ICBG742.bed', col_names = c('chr_q', 'start_q', 'end_q', 'gene_q'))

block_bed_r <- left_join(ref_bed, mcscan_blocks, by = 'gene')
block_bed_q <- left_join(q_bed, mcscan_blocks, by = 'gene_q')


