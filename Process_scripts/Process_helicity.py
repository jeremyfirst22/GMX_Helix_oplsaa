#!/usr/bin/env python

projectDir='/Users/jeremyfirst/GMX_Helix_oplsaa'
#rcFile='poster.rc'
rcFile='paper.rc'
inFile='dssp/helen.nrt'
saveDir='figures'
outname='%s/combined_helicity.pdf'%saveDir

figCols=2
figRows=3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import glob
import os
import sys
from matplotlib import rc_file

def Usage():
    print "Usage: %s <binSize=10>"%(sys.argv[0])

binSize = 10
try :
    binSize = int(sys.argv[1])
except :
    Usage()
    print "No bin size indicated, defaulting to %i"%binSize

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir) 

fig, axarr = plt.subplots(figRows, figCols, sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.1) 
fig.subplots_adjust(hspace=0.25) 
fig.text(0.5,0.05, "Time (ns)", ha='center', va='center') 
fig.text(0.05,0.5, r"Helical fraction",ha='center',va='center',rotation='vertical') 
fig.subplots_adjust(left=0.1, bottom=0.1,right=0.80,top=0.9)

blue = mpatches.Patch([],[],color='b',label=r"$\alpha$-helix")
green= mpatches.Patch([],[],color='g',label=r"$3_{10}$-helix")
red  = mpatches.Patch([],[],color='r',label=r"unfolded")

xmax=50
for row,solvent in enumerate(['water','tert','sam']): 
    for col,state in enumerate(['folded','unfolded']) : 
        indFig = plt.figure() 
        indAx = indFig.add_axes([0.1, 0.1, 0.7, 0.8])
        indFig.text(0.5,0.05, "Time (ns)", ha='center', va='center') 
        indFig.text(0.05,0.5, r"Helical fraction",ha='center',va='center',rotation='vertical') 
        ax = axarr[row,col]
        datafile = '%s_%s/%s'%(state,solvent,inFile) 
        try : 
            data = np.genfromtxt(datafile) 
        except IOError : 
            print "No file found for %s %s"%(state,solvent) 
            continue 

        frameNumber = len(data)
        assert frameNumber > 0 

        binNumber = frameNumber/binSize

        y1 = data[:,1] # Alpha 
        y2 = data[:,2] # 3_10 Helix
        y3 = data[:,3] ##unfolded
        
        y1 /= 18 
        y2 /= 18 
        y3 /= 18 

        while len(y1) %binSize !=0 : 
            y1 = y1[:-1] 
            y2 = y2[:-1] 
            y3 = y3[:-1] 
        assert len(y1) == len(y2) == len(y3) 

        y1binned = np.mean(y1.reshape(-1,binSize),axis=1) 
        y2binned = np.mean(y2.reshape(-1,binSize),axis=1) 
        y3binned = np.mean(y3.reshape(-1,binSize),axis=1) 
        
        #rc_file('%s/rc_files/%s'%(projectDir,rcFile) ) 
        
        x = data[:,0] 
        x = np.linspace(0,np.max(data[:,0]),binNumber) 
        x /= 1000 
        if x.max() > xmax : 
            xmax = x.max()
        
        ax.fill_between(x,0                ,y3binned         ,facecolor='r',linewidth=0.0)
        ax.fill_between(x,y3binned         ,y3binned+y2binned,facecolor='g',linewidth=0.0) 
        ax.fill_between(x,y3binned+y2binned,1                ,facecolor='b',linewidth=0.0) 

        indAx.fill_between(x,0                ,y3binned         ,facecolor='r',linewidth=0.0)
        indAx.fill_between(x,y3binned         ,y3binned+y2binned,facecolor='g',linewidth=0.0) 
        indAx.fill_between(x,y3binned+y2binned,1                ,facecolor='b',linewidth=0.0) 
        
        ax.set_xlim(0,np.max(x)) 
        ax.set_ylim(0,1)
        ax.set_title("%s %s"%(state,solvent)) 

        indAx.set_xlim(0,xmax) 
        indAx.set_ylim(0,1)
        indAx.set_title("%s %s"%(state,solvent)) 
    
        indFig.legend( [blue,green,red],[r"$\alpha$-helix",r"$3_{10}$-helix",r"unfolded"], 
            loc = 'center', bbox_to_anchor=(0.90, 0.5),
            fontsize='medium') 
        indFig.savefig("%s/helicity_%s_%s.pdf"%(saveDir,state,solvent), format='pdf')
        plt.close(indFig) 
        print "%9s %9s plot complete"%(state,solvent)  

#fig.legend(handles=[blue, green, red],loc=4)
fig.legend( [blue,green,red],[r"$\alpha$-helix",r"$3_{10}$-helix",r"unfolded"], 
    loc = 'center', bbox_to_anchor=(0.90, 0.5),
    fontsize='medium') 

fig.savefig(outname, format='pdf')
plt.close()

print "Combined plot complete." 




