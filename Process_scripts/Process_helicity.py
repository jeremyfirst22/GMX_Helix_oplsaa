#!/usr/bin/env python

projectDir='/Users/jeremyfirst/GMX_Helix_oplsaa'
#rcFile='poster.rc'
rcFile='paper.rc'
inFile='dssp/helen.nrt'
saveDir='figures'

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import glob
import os
import sys
from matplotlib import rc_file

def Usage():
    print "Usage: %s <Data directory> <binSize>"%(sys.argv[0])

binSize = 1000
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

datafiles = glob.glob("%s/%s"%(molec,inFile) ) 
curdir = os.path.abspath(os.getcwd())

outname = "%s/binned_%s.pdf"%(saveDir,molec) 
outname = os.path.join(curdir,outname)

for filename in datafiles: 
    print "Reading %s..."%filename 
    datai = np.genfromtxt(filename,skip_header=29,skip_footer=1)
    try : data = np.concatenate([data,datai])
    except NameError : data = datai
    print "\t Sizes: %i (%i)" %(len(datai), len(data))

frameNumber = len(data)
assert frameNumber > 0 

binNumber = frameNumber/binSize

print "Done reading data...binning" 

y1 = data[:,1] # Alpha 
y2 = data[:,2] # 3_10 Helix
y3 = data[:,3] ##unfolded

y1 /= 18 
y2 /= 18 
y3 /= 18 

y1binned =[] 
y2binned =[] 
y3binned =[] 


print frameNumber, " frames"
print binSize, " frames per bin"
print binNumber, "bins to be graphed"

for j in range(binNumber):
    sum1 = 0
    sum2 = 0
    sum3 = 0

    for i in range(binSize):
        sum1 += y1[i + j*binSize]
        sum2 += y2[i + j*binSize] 
        sum3 += y3[i + j*binSize] 

    y1binned.append(sum1/binSize)
    y2binned.append(sum2/binSize)
    y3binned.append(sum3/binSize)

y1binned = np.array(y1binned)
y2binned = np.array(y2binned)
y3binned = np.array(y3binned)

print "Data binned.... now plotting"
#rc_file('%s/rc_files/%s'%(projectDir,rcFile) ) 

x = data[:,0] 
x = np.linspace(0,np.max(data[:,0]),binNumber) 
x /= 1000 
#x = np.linspace(0,binSize*binNumber/ 10,binNumber) 

plt.fill_between(x,0                ,y3binned         ,facecolor='r',linewidth=0.0)
plt.fill_between(x,y3binned         ,y3binned+y2binned,facecolor='g',linewidth=0.0) 
plt.fill_between(x,y3binned+y2binned,1                ,facecolor='b',linewidth=0.0) 

print "Data plotted.... now labeling"

plt.xlim(0,np.max(x)) 
plt.ylim(0,1)

plt.xlabel(r"time (ns)")
plt.ylabel(r"\% helical structure")

blue = mpatches.Patch([],[],color='b',label=r"$\alpha$-helix")
green= mpatches.Patch([],[],color='g',label=r"$3_{10}$-helix")
red  = mpatches.Patch([],[],color='r',label=r"unfolded")
plt.legend(handles=[blue, green, red],loc=4)

print "Plot labeled.... now saving"

plt.savefig(outname, format='pdf')
plt.close()
print "plot saved to ", outname





