#purpose: take the most mutated intron from light chain clade of interest, isolate members that have the correct intronic divergence (based on a given number of mutations) that share
#the most mutations with this sequence as the light chain
#usage: python pick_light_chain.py clade_fasta_file mutated_light_intron_fasta germ_intron_fasta target_mutation_load (integer)

#####################################################################################
#IMPORT MODULES
#####################################################################################
import sys
import collections
from collections import defaultdict
import Bio
from Bio import SeqIO
from Bio import AlignIO
from subprocess import call
import random
import operator
import os
from Bio.Seq import Seq
import subprocess
from subprocess import call
from Bio import Phylo
from decimal import Decimal
import itertools
from itertools import *
#####################################################################################
#STORE FASTA
#####################################################################################

sequence_database={}

antibody_database={}

germline_database={}

antibody_order={}

def extract_fasta(file, dictionary):

	for seq_record in SeqIO.parse(file, "fasta"):
                dictionary[seq_record.description]=seq_record.seq.upper()

extract_fasta(sys.argv[1], sequence_database)
#store the sequence of the antibody of interest
extract_fasta(sys.argv[2], antibody_database)
extract_fasta(sys.argv[3], germline_database)
####################################################################################
#TAKE AB OF INTEREST INTRON AND ALIGN TO GERMLINE AND IDENTIFY POINT MUTATIONS
####################################################################################

target_mutation_load=int(sys.argv[4])

point_mutations=[]
deletions=[]
insertions=[]

def parse_antibody_of_interest(inputfile, point_mutations, deletions, insertions):

	num_list=['0','1','2','3','4','5','6','7','8','9']
        base_list=['A','G','C','T']

	for line in inputfile:
		stripline=line.rstrip()
		split_line=line.split("\t")
		position=int(split_line[1])
		mutation=split_line[4]
		
		i=0
		if mutation.find("-")!=-1:
			for char in mutation:
				if char=='-' and mutation[i-1]!='^':

					number=[]
					clip=mutation[(i+1):]
					num=clip[:2]
					for item in num:
						if item in num_list:
							number.append(item)
					number=int(''.join(number))
					deletions.append(position)
					deletions.append(position+number)
			
		elif mutation.find("+")!=-1:
			for char in mutation:
				if char=='-' and mutation[i-1]!='^':

                                	number=[]
                                	clip=mutation[(i+1):]
                                	num=clip[:2]
                                	for item in num:
                                        	if item in num_list:
                                        	        number.append(item)
                                	number=int(''.join(number))
                                	insertions.append(position)
                                	insertions.append(position+number)
		elif mutation.find("-")==-1 and mutation.find("+")==-1:
			for char in mutation:
				if char in base_list:
					point_mutations.append(str(position)+"-"+char)


	print point_mutations

def bowtie2_align(point_mutations,deletions,insertions, j):

	#make bowtie2 index for antibody intron of interest

	if j==2:
	#make bowtie2 index for germline intron
		call(["bowtie2-build", "{0}".format(sys.argv[3]), "germline_intron"])


	#align antibody intron to germline
	call(["bowtie2", "-f", "-x", "germline_intron", "-U", "{0}".format(sys.argv[j]), "--very-sensitive", "--end-to-end", "-S", "{0}_sam_output.txt".format(sys.argv[j]), "--un", "{0}_unmapped_reads.txt".format(sys.argv[j])])
	
	#convert the samfile to bamfile, sort it and then convert it to pileup format which is more easily parsed and comprehensible (each column is for a sequences against the reference, a dot is a match
	#a comma is a match on the reverse strand, a capital ACGT is a mismatch on the forward strand, lower case acgt is a mismatch on the reverse strand, +[0-9]+[ACGTNacgtn]+ is an insertion, -[0-9]+[ACGTNacgtn] is a deletion
	#where letters and numbers indicate size and sequence of the indel in that sequence 

	call(["samtools", "view","-b", "-o","{0}_sam_output.bam".format(sys.argv[j]), "{0}_sam_output.txt".format(sys.argv[j])])
        call(["samtools", "sort", "-T", "temp", "-o", "{0}_sam_output.sorted.bam".format(sys.argv[j]), "{0}_sam_output.bam".format(sys.argv[j])])

        subprocess.call(["samtools", "mpileup", "-f", "germline_intron.fasta","-d", "10000", "-o", "{0}_pileup.txt".format(sys.argv[j]), "{0}_sam_output.sorted.bam".format(sys.argv[j])])

	inputfile=open("{0}_pileup.txt".format(sys.argv[j]))
	
	#PARSE THE PILEUP FILE TO GET THE INTRONIC MUTATIONS OF INTEREST
	if j==2:
		parse_antibody_of_interest(inputfile, point_mutations, deletions, insertions)
	inputfile.close()	
bowtie2_align(point_mutations,deletions,insertions, 2)

	
#####################################################################################
#
#####################################################################################

#prep introns of all the sequences for alignment and print into a file

def prep_introns(sequence_database):

	outputfile_intron=open("{0}_prepped_introns.fasta".format(sys.argv[1]),"w")

	for key in sequence_database:
		split_key=key.split("-")
		antibody_length=int(split_key[2])
		seq=sequence_database[key]
		intron=seq[antibody_length:]
		outputfile_intron.write(">{0}\n".format(key))
		outputfile_intron.write("{0}\n".format(intron))
	
	outputfile_intron.close()
prep_introns(sequence_database)


######################################################################################
#
######################################################################################
#final_list=[]
def parse_abs(inputfile, point_mutations, germline_database, target_mutation_load):

	final_list=[]

	
	num_list=['0','1','2','3','4','5','6','7','8','9']
        base_list=['A','G','C','T']


	for key in germline_database:
		germline=germline_database[key]




	max_coverage=0

	for line in inputfile:
                count_mutations=0
                count_deletions=0
                count_insertions=0

     		if line.find("@")==-1 and line.find("MD")!=-1:

                        split_line=line.split("\t")
                        sequence_name=split_line[0]
                	for item in split_line:
                              	if item.find("MD:")!=-1:
                                        MD_field=item
                                        MD_split=MD_field.split(":")
                                        mutations=MD_split[2]
                                       	for item in mutations:
                                               	if item in base_list:
                                                        count_mutations+=1
                                                        #THIS INCLUDES MUTATIONS AND DELETED BASES

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


			sequence=split_line[9] #store the sequence itself
			sum=count_mutations+count_insertions+count_deletions	

			####USE THE CIGAR STRING FROM THIS SAME FILE, IT IS MUCH EASIER THAN TRYING TO USE THE PILEUP FORMAT AGAIN
			k=0			
			if sum==target_mutation_load or sum==(target_mutation_load-1):

				#begin here with new algorithm 
				for item in point_mutations:

					loc=int(item.split('-')[0])
					mut=item.split('-')[1]

					temp_num=[]
					
					final_num=[]
					final_sections=[]
					
					for char in cigar:
						if char in num_list:
							temp_num.append(char)
						if char=='M':
							match_mismatch=int(''.join(temp_num))
							final_num.append(match_mismatch)
							final_sections.append('M')
							temp_num=[]
						if char=='I':
							insertion=int(''.join(temp_num))
							final_num.append(insertion)
							final_sections.append('I')
							temp_num=[]
						if char=='D':
							deletion=int(''.join(temp_num))
							final_num.append(deletion)
							final_sections.append('D')
							temp_num=[]
							
					i=0	
					tally=0
                                        found=False
                                        for item in final_sections:
                                                if item=='M':
                                                        tally=tally+final_num[i]
                                                        if i==0 and loc<=tally: #this is if the mutation is satisifed in the first part of the sequence
                                                                if sequence[loc-1]==mut:
									k+=1

                                                                        found=True
                                                        if found==False and i>0 and loc<=tally: #check to see if in another section
                                                                total=0
                                                                for j in range(i):
                                                                        if final_sections[j]=='M':
                                                                                total=total+final_num[j]
                                                                        if final_sections[j]=='D':
                                                                                total=total-final_num[j]
                                                                        if final_sections[j]=='I':
                                                                                total=total+final_num[j]

								for w in range(i):

                                                                        if final_sections[w]=='D':
                                                                                loc=loc-final_num[w]
                                                                        if final_sections[w]=='I':
                                                                                loc=loc+final_num[w]

                                                                if len(sequence)>loc:
                                                                        if loc>=total and sequence[loc-1]==mut:
                                                                                k+=1

                                                                                found=True
                                                if item=='D':
                                                        tally=tally-final_num[i]
                                                if item=='I':
                                                        tally=tally+final_num[i]
                                                i+=1

				
			if k>0 and k==max_coverage:
				if sequence_name not in final_list: 
					final_list.append(sequence_name)
			if k>max_coverage:
				if sequence_name not in final_list:
					max_coverage=k
					final_list=[]
					final_list.append(sequence_name)



			print final_list

			global final_list
			
def bowtie2_align_repertoire(point_mutations):

	

	call(["bowtie2", "-f", "-x", "germline_intron", "-U", "{0}_prepped_introns.fasta".format(sys.argv[1]), "--end-to-end", "-S", "{0}_prepped_introns.fasta_sam_output.txt".format(sys.argv[1]), "--un", "{0}_prepped_introns.fasta_unmapped_reads.txt".format(sys.argv[1])])

 
	inputfile=open("{0}_prepped_introns.fasta_sam_output.txt".format(sys.argv[1]))	
	#parse the pileup file and look for sequences of interest

	parse_abs(inputfile, point_mutations, germline_database, target_mutation_load)

bowtie2_align_repertoire(point_mutations)
###########################################################################################
#PRINT OUT 
###########################################################################################


def print_members(final_list, sequence_database):

	outputfile=open("{0}_potential_minimal_pathway_members.fasta".format(sys.argv[1]),"w")

	if len(final_list)>1:
		random_num=random.randint(0,len(final_list))

		item=final_list[random_num]
	else:
		item=final_list[0]
	sequence=sequence_database[item]
	outputfile.write(">{0}\n".format(item))
	outputfile.write("{0}\n".format(sequence))

	outputfile.close()

print_members(final_list, sequence_database)


