For processing dataset:

Complete VRC26 light chain lineage

Scripts used in the following order:

1. productive_CDRs_FWRs_IGLV1-51_IGLJ1.py
 - retain productive sequences with 3 CDRs, 4 FWRs, annotate to IGLV1-51 and IGLJ1

2. filter_blastn_output_pident_qcovs_make_fasta.py
 - nucleotide blast for sequences with >92% nucleotide CDRL3 identity to previously published VRC26 members

3. rename_antibody_length_VL.py
 - renamve sequences to include variable region length

4. primer_matching_fixed.py
 - retain sequences with exact matches to primers

5. sort_by_prevalence.py
 - sort sequences by prevalence

6. USEARCH_cluster_centroids.py
 - cluster sequences to 96% nucleotide identity

7. make_fasta_linear.py
 - make one-liner fasta file

8. USEARCH_centroids_more_than_two_members.py
 - retain centroids from clusters >2 members

9a. VJintron_alignment.py
 - retain sequences with CDRL3 motif, translate/align/reverse translate exon segment and concatenate to intron nucleotide alignment

9b. reverse_translate_alignment.py 
 - auxiliary script for VJintron_alignment.py 