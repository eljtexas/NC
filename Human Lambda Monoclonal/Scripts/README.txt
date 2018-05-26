For processing dataset:

Human lambda light chain monoclonal control

Scripts used in the following order:

1. productive_CDRs_FWRs_IGLV1-51_IGLJ1.py
 - retain productive sequences with 3 CDRs, 4 FWRs, annotate to IGLV1-51 and IGLJ1

2. rename_antibody_length_VL.py
 - rename sequences to include variable region length

3. primer_matching_fixed.py
 - remove primers that do not have exact matches to one of the primers

4. plot_monoclonal_error_distribution.py
 - generates plots of number of mutations per sequence in sample
 - generates plots for percent of sequence with a certain number of mutations

5. sort_by_prevalence.py
 - sort sequences by prevalence

6. USEARCH_cluster_centroids.py
 - cluster sequences to 96% nucleotide identity

7. make_fasta_linear.py
 - make one-liner fasta file 

8. USEARCH_cluster_centroids.py
 - cluster sequences to 96% nucleotide identity 

9. USEARCH_centroids_more_than_two_members.py
 - retain centroids from clusters >2 members
