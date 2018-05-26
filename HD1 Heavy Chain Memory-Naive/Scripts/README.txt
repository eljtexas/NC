For processing dataset:

Healthy donor multiplex repertoire datasets

Order of script usage:

1. productive_CDRs_FWRs.py
 - retain productive sequences with 3 CDRs, 4 FWRs

2. rename_antibody_length.py
 - rename sequences to include variable region length

3. remove_sequences_with_n.py
 - remove any sequences with uncalled bases

4. primer_matching_fixed.py
 - remove sequences that do not have exact matches to one of the primers

5. sort_by_prevalence.py
 - sort sequences by prevalence

6. USEARCH_cluster_centroids.py
 - cluster sequences to 96% nucleotide identity

7. make_fasta_linear.py
 - make one-liner fasta file

8. USEARCH_centroids_more_than_two_members.py
 - retain centroids from clusters of >2 members

9. 2D_histogram_divergence_plot.py
 - generate 2D histogram of V gene and intron divergence
 - generates input files for boxplots

