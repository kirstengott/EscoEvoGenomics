for i in `ls ../../proteins/`; do diamond blastp -d db/allgenes.fa.dmnd -q ../../proteins/$i -o $i; done
for i in `ls`; do ~/scripts/parseBlast.R $i ../parsed/$i; done
