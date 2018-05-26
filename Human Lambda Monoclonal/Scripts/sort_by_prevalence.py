#purpose: take a fasta file and sort the file by prevalence of exact matches prior to clustering
#usage: python sort_by_prevalence.py fasta_file

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
#
######################################################################################

freq_database={}

def extract_fasta_freq(x, dictionary):

	for seq_record in SeqIO.parse(sys.argv[x], "fasta"):
		try:
               		dictionary[str(seq_record.seq.upper())]+=1
		except KeyError:
			dictionary[str(seq_record.seq.upper())]=1

extract_fasta_freq(1, freq_database)


sequence_database=defaultdict(list)

def extract_fasta(x, dictionary):

        for seq_record in SeqIO.parse(sys.argv[x], "fasta"):
        	dictionary[str(seq_record.seq.upper())].append(seq_record.description)
			
extract_fasta(1, sequence_database)


sorted_dict=sorted(freq_database.items(), key=operator.itemgetter(1), reverse=True)


outputfile=open("{0}_sorted_prevalence.fasta".format(sys.argv[1]),"w")

for item in sorted_dict:
	seq=list(item)[0]
	freq=list(item)[1]
	list_names=sequence_database[seq]
	for thing in list_names:
		outputfile.write(">{0};size={1}\n".format(thing, freq))
		outputfile.write("{0}\n".format(seq))	


	
