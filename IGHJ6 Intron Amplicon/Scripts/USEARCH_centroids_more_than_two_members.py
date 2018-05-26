#purpose: extract centroids from only clusters with two or more members from USEARCH output file
#usage: python USEARCH_centroids_more_than_two_members.py fasta_file.fasta cluster_file.uc

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



###########################################################################################
#BANK CENTROIDS FASTA FILE
###########################################################################################
sequence_database={}

def extract_fasta(x, dictionary):

        for seq_record in SeqIO.parse(sys.argv[x], "fasta"):
                dictionary[seq_record.description]=str(seq_record.seq.upper())
extract_fasta(1, sequence_database)


########################################################################################
#EXTRACT FINAL SEQUENCES FROM CLUSTERS WITH ONLY 2 OR MORE MEMBERS 
########################################################################################

size_dict={}

inputfile=open(sys.argv[2])


outputfile=open("{0}_GT2.fasta".format(sys.argv[1]),"w")

for line in inputfile:
	stripline=line.rstrip()
	split_line=stripline.split("\t")
	if split_line[0].find('C')!=-1:
		centroid_name=split_line[8]
		split_name=centroid_name.split(";")[0]
		size=int(split_line[2])
		if size>=2:
			outputfile.write(">{0}\n".format(split_name))
			outputfile.write("{0}\n".format(sequence_database[split_name]))		
		

