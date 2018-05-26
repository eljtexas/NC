#purpose: take a fasta file and remove sequences with n's
#usage: python remove_sequences_with_n.py fasta_file

#####################################################################################$
#IMPORT MODULES
#####################################################################################$

import sys
import collections
from collections import defaultdict
import Bio
from Bio import SeqIO
from Bio import AlignIO
from subprocess import call
from decimal import Decimal
import os
import subprocess
import matplotlib as mpl
mpl.use('Agg') #tells python to use non-interactive back-end instead of plotting to t$
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib.colors import LogNorm
import re
import scipy
from scipy import stats


total_number_indels=0
#################################################################################
# BANK FASTA FILE
#################################################################################

sequence_database=defaultdict(list)

def extract_fasta(x, dictionary):

       for seq_record in SeqIO.parse(sys.argv[x], "fasta"):
               dictionary[seq_record.description].append(str((seq_record.seq)).upper())

extract_fasta(1, sequence_database)

print "FASTA FILE BANKED"

#################################################################################
#################################################################################

outputfile=open("{0}_no_n.fasta".format(sys.argv[1]),'w')


for key in sequence_database:
	sequence=sequence_database[key][0]
	if sequence.find("N")==-1:
		outputfile.write(">{0}\n".format(key))
		outputfile.write("{0}\n".format(sequence))

outputfile.close()
	
