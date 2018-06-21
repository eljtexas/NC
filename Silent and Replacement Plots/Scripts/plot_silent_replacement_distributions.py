#purpose: script for plotting silent and replacement mutation distributions, indels in the CDRH3 are ignored, only sites aligned directly to a location in the UCA are included for that portion
#usage: python plot_silent_replacement_distributions.py fasta_file (where fasta file is the full lineage, aligned as AA, reverse translated and contains a UCA sequenced)


##########################################################################
#IMPORT MODULES
##########################################################################
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
from pylab import axes as ax
from matplotlib.colors import LogNorm
import pylab as P
from scipy.stats import gaussian_kde

import math

#############################################################################
#store aligned sequence fasta
#############################################################################

sequence_database={}

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
	sequence_database[seq_record.description]=str(seq_record.seq.upper())


week38={}
week48={}
week59={}
week119={}
week206={}

#set the UCA sequence
for key in sequence_database:
        if key.find("UCA")!=-1:
                UCA=sequence_database[key]
	if key[0]+key[1]=='38':
		week38[key]=sequence_database[key]
	if key[0]+key[1]=='48':
		week48[key]=sequence_database[key]
	if key[0]+key[1]=='59':
		week59[key]=sequence_database[key]
	if key[0]+key[1]+key[2]=='119':
		week119[key]=sequence_database[key]
	if key[0]+key[1]+key[2]=='206':
		week206[key]=sequence_database[key]

############################################################################
#Codon Table
###########################################################################


table={'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 
                              'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 
                              'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C', 
                              'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 
                              'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 
                              'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 
                              'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 
                              'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 
                              'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 
                              'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 
                              'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 
                              'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 
                              'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 
                              'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 
                              'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 
                              'GGG': 'G'} 
################################################################################
#loop through sequences and compare to UCA, then record and plot 
################################################################################

UCA_length=0
for char in UCA:
	if char!='-':
		UCA_length+=1


week38_silent_freq={}
week38_replacement_freq={}

week48_silent_freq={}
week48_replacement_freq={}

week59_silent_freq={}
week59_replacement_freq={}

week119_silent_freq={}
week119_replacement_freq={}

week206_silent_freq={}
week206_replacement_freq={}


def analyze(sequence_dict, silent_freq_dict, replacement_freq_dict):

	for key in sequence_dict:
		seq=sequence_dict[key]
	
		i=0
		while i<=(len(seq)-2):
			UCA_codon=UCA[i]+UCA[i+1]+UCA[i+2]
			seq_codon=seq[i]+seq[i+1]+seq[i+2]

			
			
			if UCA_codon!='---' and seq_codon!='---':

				UCA_AA=table[UCA_codon]
				seq_AA=table[seq_codon]
				#if there is a nucleotide mutation
				if UCA_codon!=seq_codon:
				#is there an amino acid change
					if UCA_AA==seq_AA:
						j=0
						for char1, char2 in zip(UCA_codon, seq_codon):
							j+=1
							if char1!=char2:

								site_on_UCA=i+j

								#make sure counting correctly for indels, must subtract out those wrong places
								count_dash=0
								k=0
								while k<=(site_on_UCA-1):
									if UCA[k]=='-':
										count_dash+=1
									k+=1
									
								site_on_UCA=site_on_UCA-count_dash				
			
								try:
									silent_freq_dict[site_on_UCA]+=1
								except KeyError:
									silent_freq_dict[site_on_UCA]=1					
			
			#now check under the condition that there is an amino acid change
					if UCA_AA!=seq_AA:
						j=0 #any mutation in a replaced codon gets counted as a replacement mutation, strict definition
						for char1, char2 in zip(UCA_codon, seq_codon):
							j+=1
							if char1!=char2:
								site_on_UCA=i+j
			
								#make sure counting correctly for indels
                                                                count_dash=0
								k=0
                                                                while k<=(site_on_UCA-1):
                                                                        if UCA[k]=='-':
                                                                        	count_dash+=1
									k+=1			

                                                                site_on_UCA=site_on_UCA-count_dash
			
			
							
								try:
									replacement_freq_dict[site_on_UCA]+=1
								except KeyError:
									replacement_freq_dict[site_on_UCA]=1
			i+=3

analyze(week38, week38_silent_freq, week38_replacement_freq)
analyze(week48, week48_silent_freq, week48_replacement_freq)
analyze(week59, week59_silent_freq, week59_replacement_freq)
analyze(week119, week119_silent_freq, week119_replacement_freq)
analyze(week206, week206_silent_freq, week206_replacement_freq)


#print week206_silent_freq
#print week206_replacement_freq

###################################################################
#Plot the graphs
###################################################################
def plot_graph(seq_dat, freq_dict, name):

        total_number_sequences=len(seq_dat) #calculate percent for the frequency based on the total number of sequences

        nuc_position=[]
        frequency=[]

        sorted_dict=collections.OrderedDict(sorted(freq_dict.items()))

        for key in sorted_dict:
                nuc_position.append(int(key))
                freq=(sorted_dict[key])/Decimal(total_number_sequences)*100
                freq=round(freq,2)
                frequency.append(freq)

	#set specific colors to specific bars
	color_list=[]

	position_list=[250,251,252,352,353,354,355,356,357,397,398,399]
	
	for item in nuc_position:
		if item not in position_list:
			color_list.append('black')
		elif item in position_list:
			color_list.append('red')




        plt.bar(nuc_position, frequency, color=color_list)

        plt.tick_params(axis='both', which='major', direction='out')
        plt.tick_params(axis='both', which='major', labelsize=20)
	plt.ylim(0,100)	
	plt.xlim(0,432)
	fig=plt.gcf()
	fig.set_size_inches(40, 5)

        fig.savefig("{0}_{1}.png".format(sys.argv[1], name))

	plt.close()
plot_graph(week38, week38_silent_freq, 'week38_silent')
plot_graph(week38, week38_replacement_freq, 'week38_replacement')

plot_graph(week48, week48_silent_freq, 'week48_silent')
plot_graph(week48, week48_replacement_freq, 'week48_replacement')

plot_graph(week59, week59_silent_freq, 'week59_silent')
plot_graph(week59, week59_replacement_freq, 'week59_replacement')

plot_graph(week119, week119_silent_freq, 'week119_silent')
plot_graph(week119, week119_replacement_freq, 'week119_replacement')

plot_graph(week206, week206_silent_freq, 'week206_silent')
plot_graph(week206, week206_replacement_freq, 'week206_replacement')

print len(week38)
print len(week48)
print len(week59)
print len(week119)
print len(week206)


#print week38_silent_freq
#print week38_replacement_freq

#print week48_silent_freq
#print week48_replacement_freq

#print week59_silent_freq
#print week59_replacement_freq

#print week119_silent_freq
#print week119_replacement_freq

print week206_silent_freq
print week206_replacement_freq


###########################################################################

