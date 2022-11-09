for i in `ls * | grep -v mix`; do cat $i | sed -e "s/component/${i%_c0.30.network}.component/"; done >pivot_table.tsv
