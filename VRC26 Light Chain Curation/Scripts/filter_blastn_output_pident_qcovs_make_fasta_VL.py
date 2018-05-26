#purpose: to take a tabulated blastn output file (with qseqid sseqid pident qcovs) and filter to 92% identity and 99% query coverage
#also take in the original fasta file and use it to match back and find the full length sequence
#usage: python filter_blastn_output_pident_qcovs_make_fasta_VL.py blastn_tab_output fasta_file_db


#################################################################################
# IMPORT MODULES
#################################################################################

import sys
import collections
from collections import defaultdict
import Bio
from Bio import SeqIO
from Bio import AlignIO
from subprocess import call
import random
import operator
from decimal import Decimal

#################################################################################
# BANK FASTA FILE
#################################################################################

sequence_database={}

def extract_fasta(x, dictionary):

       for seq_record in SeqIO.parse(sys.argv[x], "fasta"):
               dictionary[seq_record.id]=str(seq_record.seq).upper()

extract_fasta(2, sequence_database)

#################################################################################
#FILTER AND PRINT
#################################################################################

inputfile=open(sys.argv[1])

filtered_dictionary=defaultdict(list)


for line in inputfile:
	stripline=line.rstrip()
	split_line=stripline.split("\t")
	qseqid=split_line[0]
	sseqid=split_line[1]
	pident=float(split_line[2])
	qcovs=float(split_line[3])


	pident_limit=92.0
	qcovs_limit=99.0


	if pident>=pident_limit:

		if qseqid!='IGHV3-30*18':
			if qcovs>=qcovs_limit:
				if qseqid.find("VRC")!=-1:
					print qseqid
				filtered_dictionary[qseqid].append(sseqid+'\t'+str(pident)+'\t'+str(qcovs))



outputfile=open("{0}_{1}pident_{2}qcovs_corrected_make_fasta.fasta".format(sys.argv[1], pident_limit, qcovs_limit),"w")



check_dict={}

for key in filtered_dictionary:
	


	value=filtered_dictionary[key]
	sort_dict=defaultdict(list)
	for item in value:
		sort_by=item.split("\t")[0]
		sequence=sequence_database[sort_by]
		new_value=key+"\t"+item+"\t"+sequence
		sort_dict[key].append(new_value)
	sorted_dict=sorted(sort_dict.items(), key=operator.itemgetter(0))
	new_dict=dict(sorted_dict)
	for key in new_dict:

		thing=new_dict[key]
		for item in thing:
			split_item=item.split("\t")
			original_name=split_item[0]
			pident=split_item[1]
			qcovs=split_item[2]
			read_name=split_item[3]
			sequence=split_item[4]
				
			try:
				check_dict[pident]+=1
			except KeyError:
				check_dict[pident]=1
					
			
				outputfile.write(">{0}\n".format(pident))
				outputfile.write("{0}\n".format(sequence))
	

	
