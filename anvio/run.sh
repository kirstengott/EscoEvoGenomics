sudo anvi-matrix-to-newick data.txt -o tree.txt
sudo anvi-interactive -d data.txt -p profile.db --title "Orthologues" --tree tree.txt --manual --server-only
sudo anvi-import-misc-data additional-items-data.txt --target-data-table items --pan-or-profile-db profile.db
sudo anvi-import-misc-data additional-layers-data.txt --target-data-table layers --pan-or-profile-db profile.db
#anvi-export-state -p profile.db                   -s default                   -o default.json
sudo anvi-import-state -p profile.db -s default_edi.json -n default
inkscape --export-type=png --export-filename=Orthologues.png -d 300 -D -b "white" Orthologues.svg 
