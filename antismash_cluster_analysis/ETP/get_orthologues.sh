grep -p "\tCDS\t" ../../bigscape/gff3/mix_c0.30.network_component_3.gff3 | grep Note=genomeID:NQYS00000000_maker-glimmerHMM_NQYS01000001.1.cluster020.gbk,BGCid:BGC0001585.1 >cds_regions_focal_genome.gff3
bedtools getfasta -fi ../../fasta/NQYS00000000.fa -fo cds_regions_focal_genome.fa -bed cds_regions_focal_genome.gff3 -s 
blastx -query cds_regions_focal_genome.fa -db ../../gene_annotation/NQYS00000000.all.maker.proteins.fasta -evalue 0.001 -out cds_regions_focal_genome.blastx -outfmt 6
 ~/scripts/parseBlast.R cds_regions_focal_genome.blastx cds_regions_focal_genome.blastx.parsed 90 0.001
head -n 1 ../../annotation/Orthogroups.tsv >>cds_regions_focal_genome.orthogroups
for i in `cut -f 2 cds_regions_focal_genome.blastx.parsed`;  do grep $i ../../annotation/Orthogroups.tsv >>cds_regions_focal_genome.orthogroups; done
