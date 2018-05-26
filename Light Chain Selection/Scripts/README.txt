VRC26 light chain selection for reconstitution

Prerequisite steps:

1. Generate fasta file for clade of interest
2. Calculate number (integer) of intronic mutations desired to achieve correct intronic divergence from germline (%)
to match heavy and light chain divergence
3. Identify member of relevant light chain clade with most intronic mutations and furthest out in the topology

Scripts:
1. pick_light_chain.py
 - takes fasta file of members of clade of interest and most mutated intron sequence from
a member of that clade as well as the number of mutations (integer) necessary to have in the light chain
intron to have a desired intronic divergence and chooses a light chain with the most shared mutations with
the mutated clade member (ensuring selection of sequence relevant to light chain phylogeny) with the 
correct number of intronic mutations

Auxiliary files:
2. 59-005140-334_intron_upper_branch.fasta (Most mutated intron for upper clade used)
3. 59-006945-334_intron_lower_branch.fasta (Most mutated intron for lower clade used)