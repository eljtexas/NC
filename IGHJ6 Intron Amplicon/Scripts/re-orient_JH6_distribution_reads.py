#purpose: take the reads that consist of the JH6 gene and downstream JH6 intron and re-orient the reads such that the JH6 gene is on the 5' end
#usage: python re-orient_JH6_distribution_reads.py quality_filtered_fasta

#######################################################################################
#Import modules
#######################################################################################

import sys
import collections
from collections import defaultdict
import Bio
from Bio import SeqIO
from Bio import AlignIO
from subprocess import call
from decimal import Decimal
import subprocess

import matplotlib as mpl
mpl.use('Agg') #tells python to use non-interactive back-end instead of plotting to the screen
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib.colors import LogNorm
import pylab as P
######################################################################################
######################################################################################

def reverse(string):
        letters=list(string)
        letters.reverse()
        return "".join(letters)

def complement(string):
        basecomplement={'A':'T','C':'G','G':'C','T':'A','N':'N'}
        letters=list(string)
        letters=[basecomplement[base] for base in letters]
        return "".join(letters)

def reversecomplement(string):
        string=reverse(string)
        string=complement(string)
        return string

######################################################################################
#BANK FASTA FILE
######################################################################################

sequence_database={}

def extract_fasta(x, dictionary):

       for seq_record in SeqIO.parse(sys.argv[x], "fasta"):
               dictionary[seq_record.description]=seq_record.seq.upper()

extract_fasta(1, sequence_database)

######################################################################################
#
######################################################################################

def orient(sequence_database):

	#look for all sequences with the primer site and re-orient
	outfile=open("{0}_re-oriented.fasta".format(sys.argv[1]),"w")

	i=0

	for key in sequence_database:
		full_length=str(sequence_database[key])
		J_gene='ccacggtcaccgtctcctcag'.upper()
		#looking for exact matches and to re-orient the reads
		if full_length.find(J_gene)!=-1:
			i+=1

			outfile.write(">1-{0}\n".format(i))
			outfile.write("{0}\n".format(full_length))
		if reversecomplement(full_length).find(J_gene)!=-1:
			i+=1
			full_length=reversecomplement(full_length)
			outfile.write(">1-{0}\n".format(i))
			outfile.write("{0}\n".format(full_length))

	outfile.close()
			
orient(sequence_database)
