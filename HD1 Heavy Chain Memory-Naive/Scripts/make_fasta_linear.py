#purpose: take a fasta file and make it linear in sequence 
#usage: python make_fasta_linear.py fasta_file


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

sequence_database={}

def extract_fasta(x, dictionary):

       for seq_record in SeqIO.parse(sys.argv[x], "fasta"):
               dictionary[seq_record.description]=seq_record.seq.upper()

extract_fasta(1, sequence_database)

######################################################################################
#PRINT OUT FINAL FILE 
######################################################################################

def print_clustered_file(cluster_dictionary):
	
	outputfile=open("{0}_make_linear.fasta".format(sys.argv[1]), "w")	
	

	for key in sequence_database:
		split_key=key.split(";")[0]
		sequence=sequence_database[key]
		outputfile.write(">{0}\n".format(split_key))
		outputfile.write("{0}\n".format(sequence))

print_clustered_file(cluster_dictionary)
	
	
