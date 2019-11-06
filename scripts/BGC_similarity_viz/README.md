# Run on one network file
BGC_similarity_viz/bigscape_network2subgraph.py bigscape_output/network_files/2019-10-01_11-14-53_hybrids_glocal/NRPS/NRPS_c0.30.network bigscape_input





# Run on many network files
for i in `ls bigscape_output/network_files/2019-10-01_11-14-53_hybrids_glocal/*/*network`; do echo $i; BGC_similarity_viz/bigscape_network2subgraph.py $i antismash4/bigscape_input; done


# Create BGC plots
BGC_similarity_viz/make_BGC_plot_wrapper.R 