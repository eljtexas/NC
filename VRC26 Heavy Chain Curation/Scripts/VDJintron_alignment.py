#purpose: take a VDJ/intron fasta file of concatenated variable region and intron and align it, such that the coding region is aligned as AA and the
#intron is aligned as nucleotide, with the coding region being reverse translated back to a nucleotide alignment and both are concatenated
#usage: python VDJintron_alignment.py fasta_file 

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
import re
#####################################################################################
#STORE LONGITUDINAL RENAMED FASTA FILE AND CHARACTERIZED DONOR ANTIBODIES FROM FILES
#####################################################################################

def extract_fasta(file, dictionary):

       for seq_record in SeqIO.parse(file, "fasta"):
       		dictionary[seq_record.description]=seq_record.seq.upper()

#####################################################################################
#CDRH3 MOTIF SEARCHING TO IDENTIFY THE CORRECT TRANSLATION FRAME 
#####################################################################################

def cdr3_frame(antibody):

	antibody=Seq(str(antibody))
	AA_seq_total=str(antibody.translate())
	AA_seq=str(antibody.translate())[-55:]#force motif search to be at the end of the heavy chain (last 65 AAs)
	match=re.search(r"C[A-Z]{10,40}WG.G", AA_seq)
	
	if match!=None:

		if len(antibody)%3!=0:
			antibody=antibody[:-1] #remove any trailing nucleotide from the J gene that is not being used in the amino acids or part of the last few TVSS sequence (basically splice site)
		return str(antibody)
	if match==None:
		pre=antibody
		antibody=antibody[1:]

		antibody=Seq(str(antibody))
		AA_seq=str(antibody.translate())[-55:]
		match=re.search(r"C[A-Z]{10,40}WG.G", AA_seq)
		if match!=None:
			CDR3=(match.group())[:-3]
			if len(antibody)%3!=0:
				antibody=antibody[:-1] #remove any trailing nucleotide
			return str(antibody)
		if match==None:
			antibody=antibody[1:]
			antibody=Seq(str(antibody))
			AA_seq=str(antibody.translate())[-55:]
			match=re.search(r"C[A-Z]{10,40}WG.G",AA_seq)
			if match !=None:
				CDR3=(match.group())[:-3]
				if len(antibody)%3!=0:
					antibody=antibody[:-1]		
				return str(antibody)
			if match==None:
				return ''

#####################################################################################
#ADD VRC26 CHARACTERIZED ANTIBODIES INTO SEQUENCES AND SPLIT INTO MULTIPLE FILES
#####################################################################################

def split_sequences(sequence_database, round):

	i=0
	j=1
		
	outputfile_VDJ=open("{0}_prep_for_alignment_VDJ_part_{1}_rd{2}.fasta".format(sys.argv[1],j, round),"w")
	
	outputfile_intron=open("{0}_prep_for_alignment_intron_part_{1}_rd{2}.fasta".format(sys.argv[1],j, round),"w")	

	temp_antibody_dict={}
	temp_intron_dict={}	

	for key in sequence_database:
		split_key=key.split("-")

		antibody_length=int(split_key[2])

		seq=sequence_database[key]
		antibody=seq[:antibody_length] #extract antibody

		#look for the cdr3 in all frames and when it is found, return that sequence in the correct frame
		#this calls the cdr3_frame function to find the frame, trim the trailing g (or other base) from the J gene and to keep only those that have a C..[A-Z]{1,40}WG.G motif for the CDR3, which should help weed out other poor quality sequences
		framed_antibody=cdr3_frame(antibody)
		
		if framed_antibody!='': #basically if it met all of the requirements of being in frame, and having a CDR3 motif of C...WGxG it passes the test, otherwise it likely has sequencing error and just discard it 
			intron=seq[(antibody_length-1):] #extract intron with a minus one to take into account that trailing J gene nucleotide which is basically the splice site and should be present on all sequences	
			temp_antibody_dict[key]=framed_antibody
			temp_intron_dict[key]=intron


	keys_ab=list(temp_antibody_dict.keys())
	keys_intron=list(temp_intron_dict.keys())
	random.shuffle(keys_ab)
	random.shuffle(keys_intron)
	for item in keys_ab:
		outputfile_VDJ.write(">{0}\n".format(item))
		outputfile_VDJ.write("{0}\n".format(temp_antibody_dict[item]))
	for item in keys_intron:
		outputfile_intron.write(">{0}\n".format(item))
		outputfile_intron.write("{0}\n".format(temp_intron_dict[item]))
				
	outputfile_VDJ.close()
	outputfile_intron.close()
	temp_antibody_dict={}
	temp_intron_dict={}
	j+=1
						
	global final_count
	final_count=j
	
######################################################################################
#ALIGN EACH FILE WITH MAFFT IN AMINO ACIDS AND THEN REVERSE TRANSLATE TO NUCLEOTIDES
######################################################################################

def align_sequences(final_count, round):
	
	for j in range (final_count):

		prepped_VDJ_file="{0}_prep_for_alignment_VDJ_part_{1}_rd{2}.fasta".format(sys.argv[1],j, round)	

		prepped_intron_file="{0}_prep_for_alignment_intron_part_{1}_rd{2}.fasta".format(sys.argv[1],j, round)	

		if os.path.exists(prepped_VDJ_file):
			call(["python", "reverse_translate_alignment.py", "{0}".format(prepped_VDJ_file), "-prot_outfile", "{0}_protein_VDJ_alignment.fasta".format(prepped_VDJ_file), "-nuc_outfile", "{0}_nucleotide_VDJ_alignment.fasta".format(prepped_VDJ_file), "-aligner","mafft --auto"])
	
		if os.path.exists(prepped_intron_file):
			subprocess.call("mafft --auto {0} > {1}_nucleotide_intron_alignment.fasta".format(prepped_intron_file, prepped_intron_file), shell=True)

#####################################################################################
#CONCATENATE ANTIBODY AND INTRON NUCLEOTIDE ALIGNMENTS
#####################################################################################

def concatenate(final_count, round):

	partition_file=open("Partition.txt", "w")

	for j in range(final_count):
		temp_VDJ={}
		temp_intron={}
		
		VDJ_file="{0}_prep_for_alignment_VDJ_part_{1}_rd{2}.fasta_nucleotide_VDJ_alignment.fasta".format(sys.argv[1],j, round)

                intron_file="{0}_prep_for_alignment_intron_part_{1}_rd{2}.fasta_nucleotide_intron_alignment.fasta".format(sys.argv[1],j, round)

                if os.path.exists(VDJ_file):
                        
			extract_fasta(VDJ_file, temp_VDJ) 

                if os.path.exists(intron_file):
			extract_fasta(intron_file, temp_intron)
                        
			concat_file=open("concat_part{0}_rd{1}.fasta".format(j,round), "w")

			for key in temp_VDJ:
				ab=temp_VDJ[key]
				intron=temp_intron[key]
				combined=ab+intron
				concat_file.write(">{0}\n".format(key))
				concat_file.write("{0}\n".format(combined))				

#####################################################################################
#ALL FUNCTIONS CALLED 
#####################################################################################

consolidate_dict={}
sequence_database={}
antibody_database={}

def convergence(sequence_database):

	#take in the original sequences
	extract_fasta(sys.argv[1], sequence_database)

	k=1
	split_sequences(sequence_database, k)
	align_sequences(final_count, k)
	concatenate(final_count, k)		
convergence(sequence_database)

####################################################################################
#
####################################################################################
