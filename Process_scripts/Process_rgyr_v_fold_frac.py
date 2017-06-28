#!/usr/bin/env python

rcFile='paper.rc'
inFile ='rgyr/gyrate.xvg'
inFile2='dssp/helen.nrt'
saveDir='figures'

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm 
import glob
import os
import sys
from matplotlib import rc_file

def Usage():
    print "Usage: %s <Data directory> <binSize>"%(sys.argv[0])

binSize=100
try :
    molec = sys.argv[1].strip("/") 
except :
    Usage()
    sys.exit()
try : 
    binSize = int(sys.argv[2]) 
except : 
    Usage() 
    print "No bin size indicated, defaulting to %i"%binSize

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir) 

datafile = "%s/%s"%(molec,inFile)  
datafile2= "%s/%s"%(molec,inFile2) 
curdir = os.path.abspath(os.getcwd())

outname = "%s/rgyr_v_fold_frac_%s.pdf"%(saveDir,molec) 
outname = os.path.join(curdir,outname)

data = np.genfromtxt(datafile,skip_header=25) 
data2 = np.genfromtxt(datafile2) 

x = data[:-1,1] 
y = (data2[:-1,1] + data2[:-1,2]) / 18
print len(x), len(y) 
assert len(x) == len(y) 
assert len(x) % binSize == 0
assert len(y) % binSize == 0

xs = np.mean(x.reshape(-1,binSize),axis=1) 
ys = np.mean(y.reshape(-1,binSize),axis=1) 

colors = cm.brg(np.linspace(0,1,len(ys)) ) 
plt.plot(xs,ys,color='k',alpha=0.5,zorder=1 ) 
for x, y, c in zip(xs,ys, colors) : 
    plt.scatter(x ,y,color=c,edgecolor='none',s=40,alpha=1,zorder=2) 
plt.xlabel(r"Radius of gyration (nm)") 
plt.ylabel(r"Helical fraction") 
plt.xlim(0.8, 1.8) 
plt.ylim(0,1.0) 


plt.savefig(outname, format='pdf')
plt.close()

