for i in `ls ../proteins`; do ~/src/kofamscan/kofamscan-1.2.0/exec_annotation --cpu 40 ../proteins/$i --ko-list ~/src/kofamscan/ko_list -p ~/src/kofamscan/profiles/eukaryote.hal -o $i -f detail --report-unannotated; done
