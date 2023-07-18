#!/usr/bin/env python

from pylab import *
from numpy import loadtxt
from scipy.stats import t
from matplotlib.pyplot import *
from array import array
import numpy as np
import sys, getopt, os, re

def help ():
    print   ('BlockAv v. 0.1')
    print   ('Usage: blockAv -f <input file> -n <sizeblocks> OPTIONAL -p (file "blocksAv.dat" with av for each block)')

def blockAverage(datastream, maxBlockSize, isAvFile):
 
	Nobs         = len(datastream)           # total number of observations in datastream
	minBlockSize = 1;                        # min: 1 observation/block
	if maxBlockSize == 0:
		print("If size block is zero, Number observation/4 will be used")
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
#	print (" {} {}  {}  {}\n".format(NumBlocks, maxBlockSize,minBlockSize,Nobs))
	
				#
				#  blockSize is # observations/block
				#  run them through all the possibilities
				#
	if isAvFile:
		f=open("blocksAv.dat", "w+")
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
				if isAvFile:
					f.write("%f\n" % np.mean(datastream[ibeg:iend]))
                
		blockMean[blockCtr] = np.mean(obsProp)
		blockVar[blockCtr]  = np.var(obsProp)/(Nblock - 1)
		blockCtr += 1
 
	v = np.arange(minBlockSize,maxBlockSize)

#	print ("<x> = {0:f} +/- {1:f}\n".format(blockMean[-1], np.sqrt(blockVar[-1])))
 
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

opts, extrapar = getopt.getopt(sys.argv[1:],'f:n:hp') 
nsizeblock = 0
plotp = False

if (len(opts) == 0):
    raise Exception(help())

isAvFile=False
nsizeblock=0

for o,p in opts:
    if o in ['-f']:
        filename = p
        if not os.path.isfile(''.join(filename)):
            print
            print ('ERROR')
            print ('File ',filename,' does not exists')
            print
            sys.exit()
    elif o in ['-n']:
        nsizeblock = int(p)
    elif o in ['-p']:
        isAvFile = True   
    elif o in ['-h']:
        help ()
    else:
        assert False, "unhandled option"

if nsizeblock < 0 :
    print ("Input Error: sizeblock must be a nonnegative integer")
    sys.exit()


x = np.loadtxt(filename)
mean, var, dof  = blockAverage(x,nsizeblock,isAvFile)


if mean < np.finfo(float).eps and var < np.finfo(float).eps:
    delta = 0.0
else:
    xlo, xhi = StudentTCI(mean, var, dof)
    delta= abs((xhi-xlo)/2.0)
#print StudentTCI(mean, var, dof)

print ("\nError with 95% confidence interval")
print ("95% confidence <x> = {0:f} +/- {1:f}\n".format(mean, delta))













