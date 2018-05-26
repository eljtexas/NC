For processing dataset:

IGHJ6 healthy donor memory mutation distribution

Scripts used in the following order:

1. re-orient_JH6_distribution_reads.py
 - take reads and re-orient into 5'-3'

2. remove_sequences_with_n.py
 - remove sequences with uncalled bases

3. primer_matching_fixed.py
 - retain only sequences that have exact primer matches

4. sort_by_prevalence.py
 - sort sequences by prevalence

5. USEARCH_cluster_centroids.py
 - cluster sequences to 96% nucleotide identity

6. make_fasta_linear.py
 - make one-liner fasta file

7. USEARCH_centroids_more_than_two_members.py
 - retain only centroids from clusters >2 members

8. IGHJ6_mutation_distribution.py
 - plot IGHJ6 intron mutation distribution frequency 