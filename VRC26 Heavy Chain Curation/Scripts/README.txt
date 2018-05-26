For processing dataset:

VRC26 heavy chain longitudinal sequences

Scripts used in the following order:

1. productive_CDRs_FWRs_IGHV3-30_IGHJ3.py
 - retain productive sequences with 3 CDRs, 4 FWRs, annotate to IGHV3-30 and IGHJ3

2. filter_blastn_output_pident_qcovs_make_fasta_VH.py
 - nucleotide blast for sequences with >85% nucleotide CDRH3 identity to previously published VRC26 members

3. rename_antibody_length.py
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

9a. VDJintron_alignment.py
 - retain sequences with CDRH3 motif, translate/align/reverse translate exon segment and concatenate to intron nucleotide alignment

9b. reverse_translate_alignment.py
 - auxiliary script for VDJintron_alignment.py