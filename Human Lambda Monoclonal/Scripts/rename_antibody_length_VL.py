#purpose: take a longitudinal fasta file and an IMGT nucleotide file and 
#rename the sequence to have the longitudinal timepoint and the length of the VJ region
#so that it can easily be split into VJ and intron for downstream alignment and clustering
#usage: python rename_antibody_length_VL.py longtiduinal_fasta IMGT_3_nucleotide_file prefix_timepoint_name 

##############################################################################
#IMPORT MODULES
##############################################################################
import sys
import collections
from collections import defaultdict
import Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from subprocess import call
import random
import operator
import os

##############################################################################
#STORE THE LONGITUDINAL FASTA FILE (ALL TIME POINTS)
##############################################################################
sequence_database={}

def extract_fasta(x, dictionary):

       for seq_record in SeqIO.parse(sys.argv[x], "fasta"):
               dictionary[seq_record.description]=seq_record.seq.upper()

extract_fasta(1, sequence_database)

##############################################################################
#STORE APPSOMA DATABASE QUERY OF FULL LENGTH AND ANTIBODY SEQUENCE  
##############################################################################

IMGT_nucleotide_database=defaultdict(list)

wanted_items=['V-J-REGION']

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

#               print temp_list
                else: #do the following if it is not the first line in the file
                        for i in range(len(temp_list)):
                                header=split_line[1]
                                #header=header.split("_")[0]
#                               print header
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

extract_IMGT(2, IMGT_nucleotide_database, wanted_items)

print "IMGT SUMMARY FILE BANKED"


###############################################################################
# MATCH ANTIBODY AND RENAME TO MARK THE INTRON/VARIABLE REGION BREAK POINT
###############################################################################

updated_sequence_database={}

def match(sequence_database, IMGT_nucleotide_database, updated_sequence_database):

	outputfile=open("{0}_VDJ_length.fasta".format(sys.argv[1]),"w")
	k=0
	for key in sequence_database:
		k+=1
		if len(str(k))<7:
			zero_list='000000'
			num_zeros=6-(len(str(k)))
			zeros=zero_list[:num_zeros]
	
		antibody=str(IMGT_nucleotide_database[key][1])
		length=len(antibody)
#		if length%3!=0:
#			length=length-1
#			if length%3!=0:
#				length=length-1
		new_name="{0}-".format(sys.argv[3])+zeros+str(k)+"-{0}".format(length)#concatenate VDJ length and add longitudinal timepoint
#		print new_name
#		print sequence_database[key]

		outputfile.write(">{0}\n".format(new_name))
		outputfile.write("{0}\n".format(sequence_database[key]))

match(sequence_database, IMGT_nucleotide_database, updated_sequence_database)

