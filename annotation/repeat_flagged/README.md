# Methods

Repeat regions were extracted into BED files from maker annotation files for in house genomes, or from softmasked regions of genomic fasta files. Bedtools interesect was run to find overlaps between repeat regions and genomic exons.

```bash
bedtools intersect -b $repeats_file -a $exons_gff3 -wao
```

A python script was used to parse out the genes from the intersection that had >= 20% of the exonic gene body overlapping a repeat. These genes were flagged as repeats.