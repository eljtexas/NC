#purpose: plot the number of SHM events in the V gene and the intron separately and make boxplots for VRC26 longitudinal samples and calculate p-values
#usage: python make_figure1_boxplots_healthy_donor_points_density.py naive_intron_file naive_V_gene_file memory_intron_file memory_V_gene_file
#note:files are a new number on each line representing the number of SHM events in a unique sequence for each line, intron value is in one file and V gene value in another, same sequence, same line

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

###########################################################################
#PLACE DATA IN NEW LISTS
###########################################################################


naive=[]
memory=[]

def make_lists(x,longitudinal):
	
	inputfile=open(sys.argv[x], "r")
	
	for line in inputfile:
		stripline=line.rstrip()
		SHM_events=float(stripline)
		longitudinal.append(SHM_events)
		
	inputfile.close()

make_lists(1,naive)
make_lists(2,memory)

############################################################################
#PLOT BOXPLOTS AND CALCULATE P-VALUES
############################################################################

P.figure()

data=[naive,memory]

bp=P.boxplot(data, whis=1000000) #setting whiskers to an unreasonably large number forces a plot of the whiskers to represent the min and max of the data 


for i in range(len(data)):
	y=data[i]
	y=np.array(y)
	x=np.random.normal(i+1,0.08,size=len(y))
	x=np.array(x)

#calculate point density
	xy=np.vstack([x,y])
	z=gaussian_kde(xy)(xy)

#sort points by density
	idx=z.argsort()
	x,y,z = x[idx], y[idx], z[idx]
	axes=P.gca()
	axes.set_ylim([0,50])

	P.scatter(x,y,c=z,s=15, edgecolor='')

	P.boxplot(data, whis=1000000)

	print np.average(data[i])


plt.savefig("{0}_points_density.png".format(sys.argv[1]))

from scipy.stats import mannwhitneyu
from scipy.stats import levene

print levene(naive,memory)
print mannwhitneyu(naive,memory) #must multiply this number by two to get the final value for two-tailed Mann-Whitney U test  
