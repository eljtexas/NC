#purpose: take a stitched and q filtered fasta file and use IMGT info to keep only
#the productive reads with all 3 CDRs, 4 FWRs and that use IGHV3-30 and IGHJ3
#usage: python productive_CDRs_FWRs_IGHV3-30_IGHJ3.py 1_IMGT_summary_file

########################################################################################
#IMPORT MODULES
########################################################################################

import sys
import collections
from collections import defaultdict
import Bio
from Bio import SeqIO
from Bio import AlignIO
from subprocess import call
from decimal import Decimal

import matplotlib as mpl
mpl.use('Agg') #tells python to use non-interactive back-end instead of plotting to the screen
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib.colors import LogNorm

#########################################################################################
#STORE IMGT SUMMARY FILE
#########################################################################################

IMGT_summary_database=defaultdict(list)

wanted_items=['Functionality', 'V-GENE and allele','J-GENE and allele', 'CDR-IMGT lengths', 'FR-IMGT lengths', 'Orientation','Sequence']

def extract_IMGT(x, dictionary, wanted_items):
        inputfile=open(sys.argv[x])
        temp_list=[] #make temp_list for collecting names of columns
	j=1
        for line in inputfile: #loop through every item in row and add to dictionary using column name
                strip_line=line.rstrip()
                split_line=strip_line.split("\t")

                #this piece is just for checking for the first line
                if strip_line.find("Sequence")!=-1:
                        for item in split_line:
                                temp_list.append(item)


                else: #do the following if it is not the first line in the file
                        for i in range(len(temp_list)):
				header=split_line[1]

				if temp_list[i] in wanted_items: 
                                	dictionary[header].append(temp_list[i])
                                	try:
                                        	dictionary[header].append(split_line[i])
                                	except:
                                        	dictionary[header].append('') #check to see if there is anything there, i$
    		j+=1
		if j%10000==0:
			print j
		        
        inputfile.close()

extract_IMGT(1, IMGT_summary_database, wanted_items)

##########################################################################################
#FILTER SEQUENCES AND CREATE NEW FASTA FILE
##########################################################################################

outputfile=open("{0}_productive_genes_cdrs_frs.fasta".format(sys.argv[1]),"w")

def filter(dictionary, outputfile):

	for key in dictionary:
		value=dictionary[key]
		for i,j in enumerate(value):
			if j=='Functionality':
				functionality=value[i+1]
			if j=='V-GENE and allele':
				v_gene=value[i+1]
			if j=='J-GENE and allele':
				j_gene=value[i+1]
			if j=='CDR-IMGT lengths':
				cdrs=value[i+1]
			if j=='FR-IMGT lengths':
				frs=value[i+1]
			if j=='Sequence':
				sequence=value[i+1]
		
		if functionality=='productive':
			if v_gene.find("IGHV3-30")!=-1:
				if j_gene.find("IGHJ3")!=-1:  
					if cdrs.find("X")==-1:
						if frs.find("X")==-1:
							outputfile.write(">{0}\n".format(key))
							outputfile.write("{0}\n".format(sequence.upper()))


filter(IMGT_summary_database, outputfile)

