for i in `ls protein`; do ./bin/WoLFPSort/bin/runWolfPsortSummary fungi < protein/$i >output/$i; done
for i in `ls output/*`; do base=`basename $i`; wolfpsort2table.py $i partial >parsed_output/$base; done
for i in `ls protein`; do diamond blastp -d db/pepunit.lib.dmnd -q protein/$i -o output/$i; done
for i in `ls output`; do parseBlast.R output/$i parsed_output/$i 0; done
for i in `ls protein`; do ./run_targetp.sh -f protein/$i; done
for i in `cat ../do.txt`; do ~/src/WoLFPSort/bin/runWolfPsortSummary fungi < ../proteins/$i* >wolfpsort/$i; done
