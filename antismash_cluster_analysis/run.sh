
for gff in `ls *gff3`
do
    base=${gff%.gff3}
    ## with 1k windows
    for i in `cat $gff | cut -f 9 | cut -d ";" -f 4 | grep genomeID | sed -e "s/Note=genomeID://" | sed -e "s/_.*$//" | sort | uniq`; do bedtools window -a $gff -b ../gene_annotation/$i* >output_1k/$base-$i.gff3; done
    ## for 5k windows
    for i in `cat $gff | cut -f 9 | cut -d ";" -f 4 | grep genomeID | sed -e "s/Note=genomeID://" | sed -e "s/_.*$//" | sort | uniq`; do bedtools window -w 5000 -a $gff -b ../gene_annotation/$i* >output_5k/$base-$i.gff3; done
    
    ## for 10k windows
    for i in `cat $gff | cut -f 9 | cut -d ";" -f 4 | grep genomeID | sed -e "s/Note=genomeID://" | sed -e "s/_.*$//" | sort | uniq`; do bedtools window -w 10000 -a $gff -b ../gene_annotation/$i* >output_10k/$base-$i.gff3; done
done

## remove files without card overlaps
find ./ -size 0 -print -delete
