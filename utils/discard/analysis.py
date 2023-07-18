#!/usr/bin/env python

###############################################################################
# For csf, you need to include these two lines
#
# import matplotlib
#matplotlib.use('Agg')
#
#
# The library needed to be loaded by cfs3
#
# module load apps/pyprophet/2.0.4/gcc-4.8.5+python-2.7.8+numpy-1.10.4+pandas-0.17.0+numexpr-2.4.4+scipy-0.17.0+matplotlib-1.5.1+seaborn-0.7.0+scikitlearn-0.18.1+setuptools-24.0.1
#
#
#
################################################################################
# Questa c'era nella vecchia versione ma l'ho tolta non mi ricordo perche' :(
#from __future__ import print_function
#

from pylab import *
from numpy import loadtxt
from scipy.stats import t
#from matplotlib import rc
#from matplotlib import ticker
from matplotlib.pyplot import *
from array import array
import numpy as np
#import matplotlib.pyplot as plt
import sys, getopt, os, re

def help ():
    print   ('Analysis v. 0.1')


def blockAverage(datastream, isplot=True, maxBlockSize=0,outp="out.pdf"):
#	"""This program computes the block average of a potentially correlated timeseries "x", and 
#	provides error bounds for the estimated mean <x>. 
#	As input provide a vector or timeseries "x", and the largest block size.	
#	Check out writeup in the following blog posts for more:
#	http://sachinashanbhag.blogspot.com/2013/08/block-averaging-estimating-uncertainty_14.html
#	http://sachinashanbhag.blogspot.com/2013/08/block-averaging-estimating-uncertainty.html
#	"""
 
	Nobs         = len(datastream)           # total number of observations in datastream
	minBlockSize = 1;                        # min: 1 observation/block
	if maxBlockSize == 0:
		maxBlockSize = int(Nobs/4);        # max: 4 blocs (otherwise can't calc variance)
  
	if int(Nobs/maxBlockSize) < 4:
		print ('Variance needs four blocks, the number of the size block is too high')
		print ('or the number of the observations is too low. Check it better')
		print ('Tip, reduce nsizeblock')
		sys.exit()

	NumBlocks = maxBlockSize - minBlockSize   # total number of block sizes

	blockMean = np.zeros(NumBlocks)               # mean (expect to be "nearly" constant)
	blockVar  = np.zeros(NumBlocks)               # variance associated with each blockSize
	blockCtr  = 0
	print (" {} {}  {}  {}\n".format(NumBlocks, maxBlockSize,minBlockSize,Nobs))
	
				#
				#  blockSize is # observations/block
				#  run them through all the possibilities
				#
	f  = open("final_block.dat", "w+")
	for blockSize in range(minBlockSize, maxBlockSize):

		Nblock    = int(Nobs/blockSize)               # total number of such blocks in datastream
		obsProp   = np.zeros(Nblock)                  # container for parcelling block 

		# Loop to chop datastream into blocks
		# and take average
		for i in range(1,Nblock+1):
			
			ibeg = (i-1) * blockSize
			iend =  ibeg + blockSize
			obsProp[i-1] = np.mean(datastream[ibeg:iend])
			if blockSize == maxBlockSize-1:
				f.write("%f\n" % np.mean(datastream[ibeg:iend]))
                
		blockMean[blockCtr] = np.mean(obsProp)
		blockVar[blockCtr]  = np.var(obsProp)/(Nblock - 1)
		blockCtr += 1
 
	v = np.arange(minBlockSize,maxBlockSize)

	print ("<x> = {0:f} +/- {1:f}\n".format(blockMean[-1], np.sqrt(blockVar[-1])))
  
#	if isplot:
#
#		plt.subplot(2,1,1)
#		plt.plot(v, np.sqrt(blockVar),'ro-',lw=2)
#		plt.xlabel('block size')
#		plt.ylabel('std')
#
#		plt.subplot(2,1,2)
#		plt.errorbar(v, blockMean, np.sqrt(blockVar))
#		plt.ylabel('<x>')
#		plt.xlabel('block size')
#
#		plt.tight_layout()
#		plt.show()
#
#	else:
#		plt.subplot(2,1,1)
#		plt.plot(v, np.sqrt(blockVar),'ro-',lw=2)
#		plt.xlabel('block size')
#		plt.ylabel('std')
#
#		plt.subplot(2,1,2)
#		plt.errorbar(v, blockMean, np.sqrt(blockVar))
#		plt.ylabel('<x>')
#		plt.xlabel('block size')
#
#		savefig(output,bbox_inches='tight',pad_inches=0.2)

	#return v, blockVar, blockMean
#	print Nobs, maxBlockSize,Nobs/maxBlockSize-1
	return blockMean[-1], np.sqrt(blockVar[-1]), Nobs/maxBlockSize-1



# Taken from: https://stackoverflow.com/questions/17203403/student-t-confidence-interval-in-python
#
def StudentTCI(loc, scale, df, alpha=0.95):
	return t.interval(alpha, df, loc, scale)

# loc is a location parameter or the mean value mu
# Scale in this case is variance sigma



#opts, extrapar = getopt.getopt(sys.argv[1:],'f:o:h:r:') 
# starts at the second element of argv since the first one is the script name
# extraparms are extra arguments passed after all option/keywords are assigned
# opts is a list containing the pair "option"/"value"

opts, extrapar = getopt.getopt(sys.argv[1:],'f:o:n:ph') 


output='out.pdf'
nsizeblock = 0
plotp = False

if (len(opts) == 0):
    error()

for o,p in opts:
    if o in ['-f']:
        filename = p
        if not os.path.isfile(''.join(filename)):
            print
            print ('ERROR')
            print ('File ',filename,' does not exists')
            print
            sys.exit()
    if o in ['-o']:
        output = p
    if o in ['-n']:
        nsizeblock = int(p)
    if o in ['-p']:
        plotp = True
    elif o in ['-h']:
        help ()
#    else:
#        assert False, "unhandled option"

if nsizeblock == 0 :
    print ("Input Error, you must insert block size")
    sys.exit()

print("A",filename)
x = np.loadtxt(filename)


mean, var, dof  = blockAverage(x,plotp,nsizeblock,output)


if mean < np.finfo(float).eps and var < np.finfo(float).eps:
    delta = 0.0
else:
    xlo, xhi = StudentTCI(mean, var, dof)
    delta= abs((xhi-xlo)/2.0)
#print StudentTCI(mean, var, dof)

print ("Error with 95% confidence interval")
print ("95% confidence <x> = {0:f} +/- {1:f}\n".format(mean, delta))













