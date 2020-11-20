for i in `ls ../../gene_annotation/*protein* | grep -v phr | grep -v pin | grep -v psq`; do base=`basename $i`; blastp -query chaetocin.fa -db $i -outfmt 6 -evalue 0.001 -out protein_output/$base; done
./parse_blasts.R
for i in `cat parsed_blasts.tsv | cut -f 2 | grep -v subject`; do grep -h $i ../../gene_annotation/*gff; done >all_blast_gff.gff
for i in `ls ../../bigscape/gff3/`; do bedtools intersect -wb -a all_blast_gff.gff -b ../../bigscape/gff3/$i >gff_out/$i; done
find ./ -size 0 -print -delete


for i in `ls ../../antismash/as4_nocassis/gff/`; do bedtools intersect -wb -a all_blast_gff.gff -b ../../antismash/as4_nocassis/gff/$i >gff_out/$i; done
for i in `ls ../../antismash/gff_merged/*uniq.gff3`; do base=`basename $i`; bedtools intersect -wb -a all_blast_gff.gff -b $i >gff_merge/$base; done
for i in `ls ../../antismash/as4_genome_only/gff/`; do bedtools intersect -wb -a all_blast_gff.gff -b ../../antismash/as4_genome_only/gff/$i >gff_glim/$i; done
for i in `ls ../../antismash/as4_nocassis/gff/`; do bedtools intersect -wb -a all_blast_gff.gff -b ../../antismash/as4_nocassis/gff/$i >gff_nocassis/$i; done
for i in `ls`; do cut -f 2 $i | sort | uniq >>../out_ids.txt; done 
for i in `cat out_ids.txt | cut -f 2 | grep -v subject`; do grep -h $i gff/*; done >all_blast_gff.gff
for i in `ls ../../antismash/as4_nocassis/gff/`; do bedtools intersect -wb -a all_blast_gff.gff -b ../../antismash/as4_nocassis/gff/$i | cut -f 10,11,12,12,14,15,16,17,18 >gff_nocassis/$i; done
for i in `ls ../../antismash/as4_genome_only/gff/`; do bedtools intersect -wb -a all_blast_gff.gff -b ../../antismash/as4_genome_only/gff/$i | cut -f 10,11,12,13,14,15,16,17,18 >gff_glim/$i; done
