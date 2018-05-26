#purpose: plot the number of SHM events in the V gene and the intron separately and make boxplots for VRC26 longitudinal samples and calculate p-values
#usage: python make_figure1_boxplots_VRC26.py WK38 WK48 WK59 WK119 WK206 
#note: longitudinal files are a new number on each line representing the number of SHM events in a unique sequence for each line 

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


WK38=[]
WK48=[]
WK59=[]
WK119=[]
WK206=[]


def make_lists(x,longitudinal):
	
	inputfile=open(sys.argv[x], "r")
	
	for line in inputfile:
		stripline=line.rstrip()
		SHM_events=float(stripline)
		longitudinal.append(SHM_events)
		
	inputfile.close()

make_lists(1,WK38)
make_lists(2,WK48)
make_lists(3,WK59)
make_lists(4,WK119)
make_lists(5,WK206)


############################################################################
#PLOT BOXPLOTS AND CALCULATE P-VALUES
############################################################################

P.figure()

data=[WK38,WK48,WK59,WK119,WK206]

bp=P.boxplot(data, whis=1000000) #setting whiskers to an unreasonably large number forces a plot of the whiskers to represent the min and max of the data 


for i in range(len(data)):
	y=data[i]
	y=np.array(y)
	print i
	print y
	x=np.random.normal(i+1,0.08,size=len(y))
	x=np.array(x)	
#calculate point density
        xy=np.vstack([x,y])
        z=gaussian_kde(xy)(xy)
#sort points by density
        idx=z.argsort()
        x,y,z = x[idx], y[idx], z[idx]

        P.scatter(x,y,c=z,s=15, edgecolor='')
	bp0=P.boxplot(data, whis=1000000)
	for box in bp0['boxes']:
                box.set(color='black', linewidth=3)
        for box in bp0['whiskers']:
                box.set(color='black', linewidth=3)
        for box in bp0['medians']:
                box.set(color='black', linewidth=3)
        for box in bp0['caps']:
                box.set(color='black', linewidth=3)

	
	axes=P.gca()
	axes.set_ylim([0,35])


plt.savefig("{0}_density.png".format(sys.argv[1]))

from scipy.stats import mannwhitneyu
from scipy.stats import levene

print levene(WK38,WK48)
print levene(WK48,WK59)
print levene(WK59, WK119)
print levene(WK119,WK206)
a,b=mannwhitneyu(WK38, WK48)
print a
print b
print "Week 38-48 times 2"
print b*2
print mannwhitneyu(WK48,WK59)
c,d=mannwhitneyu(WK48, WK59)
print c
print d
print "Week 48-59 times 2"
print d*2
e,f=mannwhitneyu(WK59,WK119)
print e
print f
print "Week 59-119 times 2"
print f*2 
g,h=mannwhitneyu(WK119,WK206)
print g
print h
print "Week 119-206 times 2"
print h*2
