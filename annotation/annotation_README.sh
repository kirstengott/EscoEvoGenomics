for i in `ls db`; do rgi main -i db/${i} -o out/${i} --input_type protein -a DIAMOND -n 1 --include_loose --exclude_nudge --clean; done
wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V8.txt
wget http://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.07312018.fam-activities.txt
wget http://bcb.unl.edu/dbCAN2/download/Tools/hmmscan-parser.gz
export PATH="~/annotation/cazy/run_dbcan/tools:$PATH"
source activate dbcan
export PATH="/home/kgotting/annotation/cazy/run_dbcan/tools/signalp-4.1:$PATH"
python run_dbcan.py --out_dir ICBG733 ICBG733.all.maker.proteins.fasta protein
for i in `ls *protein* | sed -e "s/\..*$//"`; do cp -r $i ..; done
for i in `ls ../db`; do /dentigerumycin_nfs/kirsten_backup/my_interproscan/interproscan-5.32-71.0/interproscan.sh -i ../db/$i -f tsv --goterms --pathways; done
/dentigerumycin_nfs/kirsten_backup/my_interproscan/interproscan-5.32-71.0/interproscan.sh -i db/ICBG730.all.maker.proteins_fix.fasta -f tsv --goterms --pathways -cpu 20
