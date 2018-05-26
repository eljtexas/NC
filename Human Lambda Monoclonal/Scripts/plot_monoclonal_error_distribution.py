#purpose: take a filtered monoclonal dataset and plot the PCR/sequencing error distribution
#usage: python plot_monoclonal_error_distribution.py fasta_file

##################################################################################
#IMPORT MODULES
##################################################################################
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
import math
#import decimal


#################################################################################
#Parse Bowtie2 Output Sam File
#################################################################################


mutation_dictionary={}
insertion_dictionary={}
deletion_dictionary={}

def parse_sam(inputfile, mutation_dictionary, insertion_dictionary, deletion_dictionary):
	
	alphabet_list=["A","C","G","T"]
	global N_seq #keep track of how many sequences are being looked at 
	N_seq=0
	
	for line in inputfile:
		count_mutations=0
		count_deletions=0
		count_insertions=0
	
		if line.find("@")==-1:
			N_seq+=1
			split_line=line.split("\t")
			sequence_name=split_line[0]
			for item in split_line:
				if item.find("MD:")!=-1:
					MD_field=item
					MD_split=MD_field.split(":")
					mutations=MD_split[2]
					for item in mutations:
						if item in alphabet_list:
							count_mutations+=1 
							#THIS INCLUDES MUTATIONS AND DELETED BASES
			num_list=['1','2','3','4','5','6','7','8','9']
			cigar=split_line[5]
			char_list=defaultdict(list)
			if str(cigar).find("I")!=-1 or str(cigar).find("D")!=-1:
				i=0
				for char in cigar:
					if char=='I':
						count_insertions+=1
						#THIS COUNTS INSERTION GAP OPENINGS 
					if char=='D':
						before=cigar[:i]
						last_nums=before[-2:]
						number=[]
						for item in last_nums:
							if item in num_list:
								number.append(item)
						if len(number)>1:
							number=int(''.join(number))
							count_mutations=count_mutations-number
						elif len(number)==1:
							number=int(number[0])
							count_mutations=count_mutations-number		
						count_deletions+=1
						#THIS COUNTS DELETION GAP OPENINGS
					i+=1
				
			mutation_dictionary[sequence_name]=count_mutations
			insertion_dictionary[sequence_name]=count_insertions
			deletion_dictionary[sequence_name]=count_deletions
				


#################################################################################
#Align to Sanger Sequence with Bowtie2
#################################################################################


def bowtie_align():
	
	call(["bowtie2","-f","-x","VRC26_lambda_monoclonal","-U","{0}".format(sys.argv[1]), "--end-to-end", "--very-sensitive", "-S","{0}_sam_output.txt".format(sys.argv[1])])

	inputfile=open("{0}_sam_output.txt".format(sys.argv[1]))
		
	parse_sam(inputfile, mutation_dictionary, insertion_dictionary, deletion_dictionary)
	
	inputfile.close()

bowtie_align()	


#################################################################################
#Calculate RMS values for this dataset
#################################################################################


def calculate_RMS(mutation_dictionary, insertion_dictionary, deletion_dictionary):

	#calculate RMS_mut
	
	total_mutations=0

	for key in mutation_dictionary:
		value=int(mutation_dictionary[key])
		total_mutations=total_mutations+value
	
	mean_mutations=total_mutations/Decimal(N_seq)
	mean_mutations=round(mean_mutations,3)
	
	Inner_mutation=0

	for key in mutation_dictionary:
		value=int(mutation_dictionary[key])
		Inner_mutation=Inner_mutation+(math.pow((value-mean_mutations), 2)/(N_seq))
		
	RMS_mut=math.sqrt(Inner_mutation)

	#calculate RMS_ins
	
	total_insertions=0

        for key in insertion_dictionary:
                value=int(insertion_dictionary[key])
                total_insertions=total_insertions+value

        mean_insertions=total_insertions/Decimal(N_seq)
        mean_insertions=round(mean_insertions,3)

        Inner_insertion=0

        for key in insertion_dictionary:
                value=int(insertion_dictionary[key])
                Inner_insertion=Inner_insertion+(math.pow((value-mean_insertions), 2)/(N_seq))

        RMS_ins=math.sqrt(Inner_insertion)
	
	#calculate RMS_del
	
	total_deletions=0

        for key in deletion_dictionary:
                value=int(deletion_dictionary[key])
                total_deletions=total_deletions+value

        mean_deletions=total_deletions/Decimal(N_seq)
        mean_deletions=round(mean_deletions,3)

        Inner_deletion=0

        for key in deletion_dictionary:
                value=int(deletion_dictionary[key])
                Inner_deletion=Inner_deletion+(math.pow((value-mean_deletions), 2)/(N_seq))

        RMS_del=math.sqrt(Inner_deletion)

	length=545


        outputfile=open("{0}_RMS_distribution.txt".format(sys.argv[1]),"w")

		print "RMS_mut_normalized per 100nt"
        print round(RMS_mut/length*100, 3)
        RMS_mut_norm=round(RMS_mut/length*100,3)
        print "RMS_ins_normalized per 100nt"
        print round(RMS_ins/length*100, 3)
        RMS_ins_norm=round(RMS_ins/length*100, 3)
        print "RMS_del_normalized per 100nt"
        print round(RMS_del/length*100, 3)
        RMS_del_norm=round(RMS_del/length*100, 3)


        outputfile.write("{0}\n".format(sys.argv[1]))
        outputfile.write("RMS_mut_unnormalized\n")
        outputfile.write("{0}\n".format(RMS_mut))
        outputfile.write("RMS_ins_unnormalized\n")
        outputfile.write("{0}\n".format(RMS_ins))
        outputfile.write("RMS_del_unnormalized\n")
        outputfile.write("{0}\n".format(RMS_del))
        outputfile.write("RMS_mut_normalized per 100nt or %\n")
        outputfile.write("{0}\n".format(RMS_mut_norm))
        outputfile.write("RMS_ins_normalized per 100nt or %\n")
        outputfile.write("{0}\n".format(RMS_ins_norm))
        outputfile.write("RMS_del_normalized per 100nt or %\n")
        outputfile.write("{0}\n".format(RMS_del_norm))
        outputfile.write("Mean Number of Point Mutations Per Sequence\n")
        outputfile.write("{0}\n".format(mean_mutations))
        outputfile.write("Mean Number of Deletion Events Per Sequence\n")
        outputfile.write("{0}\n".format(mean_deletions))
        outputfile.write("Mean Number of Insertion Events Per Sequence\n")
        outputfile.write("{0}\n".format(mean_insertions))
	
	mutation_frequency_dict={}

        for key in mutation_dictionary:
                num=int(mutation_dictionary[key])
                try:
                        mutation_frequency_dict[num]+=1
                except KeyError:
                        mutation_frequency_dict[num]=1

        mutation_list=[]
        frequency_list=[]

        for key in mutation_frequency_dict:
                mutation_list.append(key)
                frequency_list.append(mutation_frequency_dict[key])
	plt.figure(1)
        plt.scatter(mutation_list,frequency_list)
        plt.title('Error Rate Distribution')
        plt.xlabel('Number of Mutations Per Sequence')
        plt.ylabel('Number of Sequences')
        plt.savefig("{0}_error_rate_distribution.png".format(sys.argv[1]))



	#sort the mutation dictionary and plot the number of percent of reads per mutation number


	import operator

	sorted_dict=sorted(mutation_frequency_dict.items(), key=operator.itemgetter(0))
	

	total_num=0
	for item in sorted_dict:
		total_num+=int(item[1])
	

	percent_list=[]
	mut_list=[]

	outputfile.write("Number of Point Mutations Per Sequence and Percent of Data with Less than that Number\n")

	num=0
	for item in sorted_dict:
		mut=item[0]
		num+=int(item[1])
		percent=num/Decimal(total_num)*100
		mut_list.append(mut)
		percent_list.append(percent)

		outputfile.write("{0}	{1}\n".format(mut,percent))

	plt.figure(2)
	plt.scatter(mut_list,percent_list)
	plt.title('Percent of Data With Less than a Certain Number of Mutations')
	plt.xlabel('Number of Mutations Per Sequence')
	plt.ylabel('Percent of Sequences in Data Set')
	plt.savefig("{0}_percent_data_error_plot.png".format(sys.argv[1]))


calculate_RMS(mutation_dictionary, insertion_dictionary, deletion_dictionary)
	
