#purpose: take amplicons of IGHJ6 introns from mBCs and plot mutation distribution across sequence
#this version utilizes mafft instead of bowtie and fixes other bugs
#usage: python IGHJ6_mutation_distribution.py fasta_file

#######################################################################################
#Import modules
#######################################################################################
import os
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

#######################################################################################
#extract information from the fasta file
#######################################################################################


sequence_database=defaultdict(list)

def extract_fasta(x, dictionary):
        i=1
        for seq_record in SeqIO.parse("{0}".format(sys.argv[1]), "fasta"):
                dictionary[i].append(seq_record.id)
                dictionary[i].append((seq_record.seq).upper())
                i+=1
extract_fasta(1, sequence_database)

#print sequence_database
print "FASTA FILE BANKED"
#######################################################################################
#use J gene from IMGT #3 information to determine where the intron starts in the fasta file sequences
#######################################################################################
record_intron={}
intron_frequency_dict={} #this was added in the collapse_introns version of this script because it was too slow, so i collapse them to unique introns and count the frequency
def isolate_introns(sequence_database, intron_frequency_dict):
        ####################################isolate all of the introns after re-orienting the sequence in the sequence database 
        total_reads=0
        count_not_found=0
        for key in sequence_database:
                full_length=str(sequence_database[key][1])
                J_gene="ccacggtcaccgtctcctcag".upper()
                if full_length.find(J_gene)!=-1:
                        start_index=full_length.find(J_gene)
                        len_J_gene=len(J_gene)
                        intron=str(full_length[(start_index+len_J_gene):])
			reverse_primer="actgaggtcctggagcctcc".upper()
			reverse_primer=reversecomplement(reverse_primer)
			#check to make sure no reverse primers outside of exact matching primer so that it doesn't affect the alignment 
			if intron.find(reverse_primer)!=-1:
				start_index=intron.find(reverse_primer)
				len_intron=len(intron)
				intron=str(intron[:(start_index+20)])
				
				#intron trimmed of any extraneous nucleotides 
                        	sequence_database[int(key)].append("Full Length")
                        	sequence_database[int(key)].append(full_length)
                        	sequence_database[int(key)].append("Intron")
                        	sequence_database[int(key)].append(intron)

                elif reversecomplement(full_length).find(J_gene)!=-1: ###################take into account both possible directions
                        full_length=reversecomplement(full_length)
                        start_index=full_length.find(J_gene)
                        len_J_gene=len(J_gene)
                        intron=str(full_length[(start_index+len_J_gene):])
			reverse_primer="actgaggtcctggagcctcc".upper()
                        reverse_primer=reversecomplement(reverse_primer)
                        #check to make sure no reverse primers outside of exact matching primer so that it doesnt affect the alignment
                        if intron.find(reverse_primer)!=-1:
                                start_index=intron.find(reverse_primer)
                                len_intron=len(intron)
                                intron=str(intron[:(start_index+20)])
				#intron trimmed of any extraneous nucleotides 
	                        sequence_database[int(key)].append("Full Length")
	                        sequence_database[int(key)].append(full_length)
	                        sequence_database[int(key)].append("Intron")
	                        sequence_database[int(key)].append(intron)
                elif reversecomplement(full_length).find(J_gene)==-1 and full_length.find(J_gene)==-1:
                        count_not_found+=1

isolate_introns(sequence_database, intron_frequency_dict)
##############################################################################################
#parse mafft file and count frequency of mutations
##############################################################################################
def parse_aln_file(inputfile, final_location_dictionary):


	alphabet_list=["A","C","T","G"]

	alignment=AlignIO.read(inputfile,"clustal")

	germ=str(alignment[0].seq)

	second=str(alignment[1].seq)
	
	location_on_germ=0
	
	number_point_mutations=0

	number_deletions=0

	number_insertions=0

	counter=0
	for char1, char2 in zip(germ, second): # loop over both aligned sequences 
		counter+=1
		index=counter-1
		if char1!='-':	
			location_on_germ+=1 #if no insertion in our sequence and no gap in germline, count as original location
			if char2!='-': #if no deletion
				if char1!=char2:
					final_location_dictionary[location_on_germ]+=1 # keep track in global dictionary at this specific site 
					number_point_mutations+=1
					
			if char2=='-':
				if (index-1)>=0:
					if second[index-1]!='-':
						final_location_dictionary[location_on_germ]+=1
						number_deletions+=1 #this counts at the deletion start site

		if char1=='-':
			if (index-1)>=0:
				if germ[index-1]!='-':
					location_on_germ+=1
					#temporarily marked the insertion as occuring on the next site
					final_location_dictionary[location_on_germ]+=1
					number_insertions+=1
					#return the location on germ to what it was previously
					location_on_germ-=1


##############################################################################################
#mafft alignment 
##############################################################################################


final_location_dictionary={}

def mafft_align(sequence_database, final_location_dictionary):

	#populate original location dictionary with all zero SHM events at every position

	germline_file=open("human_jh603_intron_dist.fasta","r")
	

	for line in germline_file:
		stripline=line.rstrip()
		if stripline.find(">")==-1:
			germ_length=len(stripline)

	for i in range(germ_length):
		i+=1
		final_location_dictionary[i]=0
		

	#prepare alignment final and run mafft 
	i=0
	for key in sequence_database:
		i+=1
		if i%100==0:
			print i
		value=sequence_database[key]
		for item in value:
			if str(item).find("Intron")!=-1:
				index=value.index(item)
				data_index=index+1
				intron=value[data_index]
				
				germline_file=open("human_jh603_intron_dist.fasta", "r")

				with open("temp_output_pairwise_mafft.fasta","w") as temp_output_pairwise_alignment:
					
					for line in germline_file:
						stripline=line.rstrip()
						temp_output_pairwise_alignment.write("{0}\n".format(stripline))

					
					temp_output_pairwise_alignment.write(">{0}\n".format(key))
					temp_output_pairwise_alignment.write("{0}\n".format(intron))
					temp_output_pairwise_alignment.close()
					subprocess.call('mafft --auto --quiet --clustalout temp_output_pairwise_mafft.fasta > output.aln', shell=True)

			
					inputfile=open("output.aln")
					parse_aln_file(inputfile, final_location_dictionary)
					inputfile.close()
					

					os.remove("temp_output_pairwise_mafft.fasta")
					germline_file.close()	
	

mafft_align(sequence_database, final_location_dictionary)

################################################################################################
#Plot intronic mutation distribution 
################################################################################################
def plot_graph(sequence_database, final_location_dictionary):

	total_number_sequences=len(sequence_database)

	nuc_position=[]
	frequency=[]
		
	sorted_dict=collections.OrderedDict(sorted(final_location_dictionary.items()))

	for key in sorted_dict:
		nuc_position.append(int(key))
		freq=(sorted_dict[key])/Decimal(total_number_sequences)*100
		freq=round(freq,2)
		frequency.append(freq)

	plt.bar(nuc_position, frequency)
	
	plt.tick_params(axis='both', which='major', direction='out')
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.savefig("{0}_mutation_distribution_mafft.png".format(sys.argv[1]))		

plot_graph(sequence_database, final_location_dictionary)

###################################################################################################
#
###################################################################################################
