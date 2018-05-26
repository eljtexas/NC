#purpose: take a fasta file and cluster it with USEARCH to 96% identity 
#usage: python USEARCH_cluster_centroids.py fasta_file


######################################################################################
#IMPORT MODULES
######################################################################################
import sys
import collections
from collections import defaultdict
import Bio
from Bio import SeqIO
from Bio import AlignIO
from subprocess import call
import random
import operator

######################################################################################
#CLUSTER USEARCH
######################################################################################
def USEARCH_cluster():
	
	inputfile="{0}".format(sys.argv[1])

	outfile="{0}_USEARCH_output_96.uc".format(sys.argv[1])
	
	call(["usearch","-cluster_smallmem","{0}".format(inputfile), "-uc","{0}".format(outfile), "-id","0.96", "-sortedby", "size", "-maxaccepts","0", "-maxrejects", "0", "-centroids", "{0}_centroids.fasta".format(outfile)])
	
USEARCH_cluster()


	
	
