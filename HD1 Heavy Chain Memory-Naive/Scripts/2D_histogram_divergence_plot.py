#purpose: take a fasta file of sequences and IMGT files (#3,#8) to plot the mutational divergence of the V gene (prior to the CDR3) against the intronic mutational divergence, where indels are assumed to be single mutations (events) (in reality the divergence may be higher in these cases)
#usage: python 2D_histogram_divergence_plot.py fasta_file 1_IMGT_file 3_IMGT_file 8_IMGT_file 


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

#########################################################################################
#STORE IMGT SUMMARY FILE
#########################################################################################

IMGT_summary_database=defaultdict(list)

wanted_items=['Sequence']

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

extract_IMGT(2, IMGT_summary_database, wanted_items)

print "IMGT SUMMARY FILE BANKED"

#########################################################################################
#STORE IMGT SUMMARY FILE
#########################################################################################

IMGT_nucleotide_database=defaultdict(list)

wanted_items=['Functionality', 'V-GENE and allele','J-GENE and allele', 'J-REGION']

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

extract_IMGT(3, IMGT_nucleotide_database, wanted_items)

print "IMGT NUCLEOTIDE FILE BANKED"

#########################################################################################
#STORE IMGT SUMMARY FILE
#########################################################################################

IMGT_mutation_database=defaultdict(list)

wanted_items=['V-REGION Nb of nucleotides','V-REGION Nb of mutations', 'V-REGION Nb of silent mutations','V-REGION Nb of nonsilent mutations']

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

extract_IMGT(4, IMGT_mutation_database, wanted_items)

print "IMGT MUTATION FILE BANKED"

###################################################################################
#ISOLATE INTRONS
###################################################################################

def isolate_introns(sequence_database,IMGT_summary_database,IMGT_nucleotide_database, IMGT_mutation_database):

	outputfile=open("{0}_prepped_introns.fasta".format(sys.argv[1]),"w")

	for key in sequence_database:
		nucleotide_value=IMGT_nucleotide_database[key]
		summary_value=IMGT_summary_database[key]		

		sequence_start=summary_value.index('Sequence')
		sequence=summary_value[sequence_start+1]

		j_gene_start=nucleotide_value.index('J-REGION')
		j_gene=nucleotide_value[j_gene_start+1]
		
		j_gene_start_sequence=sequence.index(j_gene)
		len_j_gene=len(j_gene)
		intron=str(sequence[(j_gene_start_sequence+len_j_gene):])
				

		IMGT_summary_database[key].append(intron)
		
		#prepare intron fasta file for downstream alignment with bowtie2
		outputfile.write(">{0}\n".format(key))
		outputfile.write("{0}\n".format(intron))

	
isolate_introns(sequence_database, IMGT_summary_database, IMGT_nucleotide_database, IMGT_mutation_database)		

####################################################################################
#V GENE DIVERGENCE
####################################################################################

def v_gene_divergence(sequence_database, IMGT_summary_database, IMGT_nucleotide_database, IMGT_mutation_database):

	for key in sequence_database:
		mutation_value=IMGT_mutation_database[key]
		
		total_mut_start=mutation_value.index('V-REGION Nb of mutations')
		total_mut=int(mutation_value[total_mut_start+1].split(" ")[0])
		
		total_nt_start=mutation_value.index('V-REGION Nb of nucleotides')
		total_nt=int(mutation_value[total_nt_start+1].split(" ")[0])	
	
		if total_mut!='':
			v_gene_divergence=round(total_mut/Decimal(total_nt)*100, 2)
			IMGT_summary_database[key].append(v_gene_divergence)
		else:
			IMGT_summary_database[key].append('*{0}*'.format(total_mut))


v_gene_divergence(sequence_database, IMGT_summary_database, IMGT_nucleotide_database, IMGT_mutation_database)

####################################################################################
#PARSE SAM OUTPUT FILE AND CALCULATE INTRON DIVERGENCE
####################################################################################

def parse_aln_file(inputfile, IMGT_summary_database, j_gene_allele):
	alphabet_list=["A","C","T","G"]
	
	alignment=AlignIO.read(inputfile, "clustal")
	
	germ=str(alignment[0].seq)	
	second=str(alignment[1].seq)
	
	match=re.findall(r'-*', germ)
	
	num_indel_events=0
	for item in match:
		if len(item)>0:
			num_indel_events+=1

	match=re.findall(r'-*', second)

	for item in match:
		if len(item)>0:
			num_indel_events+=1	

	#NOW COUNT POINT MUTATIONS
	num_point_mutations=0

	for char1, char2 in zip(germ, second): #loop over both aligned sequences and see if any bases do not match, if so record and calculated divergence values for the sequence 
		if char1!='-':
			if char2!='-':
				if char1!=char2:
					num_point_mutations+=1

	total_mutations=num_indel_events+num_point_mutations		

	if j_gene_allele=="IGHJ1":
        	total_nt=104
      	if j_gene_allele=="IGHJ2":
             	total_nt=121
      	if j_gene_allele=="IGHJ3":
              	total_nt=117
      	if j_gene_allele=="IGHJ4":
             	total_nt=109
      	if j_gene_allele=="IGHJ5":
             	total_nt=117
     	if j_gene_allele=="IGHJ6":
             	total_nt=114
							
	name=str(alignment[1].description)		
	intron_divergence=total_mutations/Decimal(total_nt)*100
	intron_divergence=round(intron_divergence,2)
	IMGT_summary_database[name].append(intron_divergence)
	if num_indel_events>0:
		IMGT_summary_database[name].append(1)
	elif num_indel_events==0:
		IMGT_summary_database[name].append(0)

####################################################################################
#BOWTIE2 INTRON ALIGNMENT
####################################################################################

def bowtie2align(IMGT_summary_database, total_number_indels):
	
	j_dict={"IGHJ1":"human_jh1_intron", "IGHJ2":"human_jh2_intron", "IGHJ3":"human_jh302_intron", "IGHJ4":"human_jh402_var2_intron", "IGHJ5":"human_jh502_intron","IGHJ6":"human_jh603_intron"}

	i=0
	for key in IMGT_summary_database:
		i+=1
		if i%100==0:
			print i	
		nucleotide_value=IMGT_nucleotide_database[key]
		j_gene_allele_start=nucleotide_value.index('J-GENE and allele')
		j_gene_allele=nucleotide_value[j_gene_allele_start+1]
		split_allele=j_gene_allele.split("Homsap")
		if len(split_allele)>1:
			j_gene_allele=split_allele[1]
			
		if j_gene_allele.find("IGHJ1")!=-1:
			j_gene_allele="IGHJ1"
		if j_gene_allele.find("IGHJ2")!=-1:
                        j_gene_allele="IGHJ2"
		if j_gene_allele.find("IGHJ3")!=-1:
                        j_gene_allele="IGHJ3"
		if j_gene_allele.find("IGHJ4")!=-1:
                        j_gene_allele="IGHJ4"
		if j_gene_allele.find("IGHJ5")!=-1:
                        j_gene_allele="IGHJ5"
		if j_gene_allele.find("IGHJ6")!=-1:
                        j_gene_allele="IGHJ6"
		file_to_open=j_dict[j_gene_allele]

		germline_file=open("{0}.fasta".format(file_to_open))

        	with open("temp_output_pairwise_mafft.fasta", 'w') as temp_output_pairwise_alignment:

        		for line in germline_file:
                		stripline=line.rstrip()
                		temp_output_pairwise_alignment.write("{0}\n".format(stripline))
			

			value=IMGT_summary_database[key]
			intron=value[2]
			temp_output_pairwise_alignment.write(">{0}\n".format(key))
			temp_output_pairwise_alignment.write("{0}\n".format(intron))
			temp_output_pairwise_alignment.close()			
			subprocess.call('mafft --auto --quiet --clustalout temp_output_pairwise_mafft.fasta > output.aln', shell=True)	

			inputfile=open("output.aln")
			parse_aln_file(inputfile,IMGT_summary_database,j_gene_allele)
			inputfile.close()		

			os.remove("temp_output_pairwise_mafft.fasta") # delete and start over again
			germline_file.close()
		#here take into account if the intron divergence is to too high to re-align to all introns and take the highest
		intron_divergence=IMGT_summary_database[key][4]
		new_list=[]
		indel_list=[]
		if intron_divergence>0:
			for thing in j_dict:
				file_to_open=j_dict[thing]
				germline_file=open("{0}.fasta".format(file_to_open))
				with open("temp_output_pairwise_mafft.fasta", 'w') as temp_output_pairwise_alignment:

			        	for line in germline_file:
                     				stripline=line.rstrip()
                                		temp_output_pairwise_alignment.write("{0}\n".format(stripline))


                    				#value must stay as it was before, same key
                        		intron=value[2]
                        		temp_output_pairwise_alignment.write(">{0}\n".format(key))
                        		temp_output_pairwise_alignment.write("{0}\n".format(intron))
                       			temp_output_pairwise_alignment.close()
                        		subprocess.call('mafft --auto --quiet --clustalout temp_output_pairwise_mafft.fasta > output.aln', shell=True)

                        		inputfile=open("output.aln")
                        		parse_aln_file(inputfile,IMGT_summary_database,thing)
                        		inputfile.close()

                        		os.remove("temp_output_pairwise_mafft.fasta") # delete and start over again
                        		germline_file.close()
						
			new_list.append(value[4])
			new_list.append(value[6])
			new_list.append(value[8])
			new_list.append(value[10])
			new_list.append(value[12])
			new_list.append(value[14])
			new_list.append(value[16])

			indel_list.append(value[5])
			indel_list.append(value[7])
			indel_list.append(value[9])
			indel_list.append(value[11])
			indel_list.append(value[13])
			indel_list.append(value[15])
			indel_list.append(value[17])


			index_min=new_list.index(min(new_list))
			indel_boolean=indel_list[index_min]
			total_number_indels+=indel_boolean		

			if total_number_indels%100==0:
				print "number of indels"
				print total_number_indels
				print i
			value[4]=min(new_list) #replace the fourth element with the lowest intronic divergence 


bowtie2align(IMGT_summary_database, total_number_indels)

####################################################################################
#PLOT V GENE VS INTRON DIVERGENCE
####################################################################################


def divergence_plot(sequence_database, IMGT_summary_database):
	
	intron_data=[]
	v_gene_data=[]

	for key in sequence_database:
		summary_value=IMGT_summary_database[key]
		
		intron_data.append(summary_value[4])
		v_gene_data.append(summary_value[3])
		
	#print intron and v_gene data
	output1=open("{0}_intron_data.txt".format(sys.argv[1]),'w')
	for item in intron_data:
		output1.write("{0}\n".format(item))
	
	output2=open("{0}_V_gene_data.txt".format(sys.argv[1]),'w')
	for item in v_gene_data:
		output2.write("{0}\n".format(item))

	a_ma = np.ma.masked_where(intron_data > 0, intron_data)
        b_ma = np.ma.masked_where(v_gene_data < 0, v_gene_data)

        hist, xbins, ybins = np.histogram2d(a_ma,b_ma, bins=[50,50], normed=True, range=[(0,50),(0,50)])
        extent = [0,50,0,50]

        im = plt.imshow(np.ma.masked_where(hist == 0, hist).T, interpolation='none', origin='lower', extent=extent)

        cbar=plt.colorbar(im)

	x=np.arange(0,51,0.001) #makes arrays
        y=np.arange(0,51,0.001)
	xticks=np.arange(0, 51, 10)
	yticks=np.arange(0, 51, 10)
	plt.yticks(yticks)
	plt.xticks(xticks)
        plt.tick_params(axis='both', which='major', direction='out', labelsize=30)
	plt.plot(x,y, 'r', label='Linear Reference') 
	plt.savefig("test_case_indels.png")

	print total_number_indels
	print len(IMGT_summary_database)

divergence_plot(sequence_database, IMGT_summary_database)
			

			


