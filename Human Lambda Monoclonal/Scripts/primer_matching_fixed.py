#purpose: take a fast file and filter for only sequences that contain exact matches to the primer sequences for both forward and reverse to help increase quality of reads
#usage: python primer_matching.py fwd_primers_fasta rev_primers_fasta fasta_to_filter

######################################################################################
#IMPORT MODULES
######################################################################################

import sys
import collections
from collections import defaultdict
import Bio
from Bio import SeqIO
from Bio import AlignIO
from subprocess import call
import random
import operator

######################################################################################
#STORE FASTA FILE
######################################################################################


forward_dict={}
reverse_dict={}

sequence_database={}

def extract_fasta(x, dictionary):

        for seq_record in SeqIO.parse(sys.argv[x], "fasta"):
        	dictionary[seq_record.description]=str(seq_record.seq.upper())

extract_fasta(1, forward_dict)
extract_fasta(2, reverse_dict)
extract_fasta(3, sequence_database)

#######################################################################################
#MATCH AND FILTER
#######################################################################################


outputfile=open("{0}_primer_match_filtered_fixed.fasta".format(sys.argv[3]),"w")

for key in sequence_database:
	sequence=sequence_database[key]
	i=0
	for key1 in reverse_dict: #loop through intron primers (should just be one)
		r_primer=reverse_dict[key1]
		if r_primer in sequence: #if you find the intron primer, keep going 
			i+=1
			break
	for key2 in forward_dict: #loop through FR1 primers
		if i==1:		#if none has been found yet keep going
			f_primer=forward_dict[key2]
			if f_primer in sequence:
				i+=1  #write to file if you found it
				break		
	if i==2:
		outputfile.write(">{0}\n".format(key))
		outputfile.write("{0}\n".format(sequence))


outputfile.close()
