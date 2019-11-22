#!/usr/bin/bash

for i in `ls *fastq.gz | cut -d _ -f 1 | uniq`;
do echo $i
   ls $i*fastq.gz | parallel gunzip {}
   jellyfish count -C -m 31 -s 5000000000 -t 16 *.fastq -o reads.jf;
   kat hist reads.jf -t 16 -o /home/kgotting/jellyfish/kat_hist_mito/$i.histo
   rm reads.jf
   ls $i*fastq | parallel gzip {}
done
