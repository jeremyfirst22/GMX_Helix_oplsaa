#!/usr/bin/env python

projectDir='/Users/jeremyfirst/GMX_Helix_oplsaa'
#rcFile='poster.rc'
rcFile='presentation.rc'
inFile='dssp/helen.nrt'
saveDir='figures'
outname='%s/bounds_helicity.png'%saveDir
outname2='%s/bounds_average_helicity.png'%saveDir
outnameD='%s/bounds_firstD_helicity.png'%saveDir

figCols=2
figRows=3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import glob
import os
import sys
from matplotlib import rc_file
import matplotlib as mpl

rc_file('rc_files/%s'%rcFile) 

mpl.rcParams['figure.figsize'] = 5,3.8

def Usage():
    print "Usage: %s <binSize=10> <endTime=-1 (ns)>"%(sys.argv[0])

binSize = 10
try :
    binSize = int(sys.argv[1])
except :
    Usage()
    print "No bin size indicated, defaulting to %i"%binSize

endTime = -1 
try : 
    endTime = int(sys.argv[2]) 
    endTime *= 1000 ##ns -> ps 
    endTime /= 4 ##ps -> frames
except : 
    print "No endTime indicated, using all available frames" 

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir) 

fig, axarr = plt.subplots(figRows, figCols, sharex='none',sharey='row', figsize=(10.5,6.5)) 
fig.subplots_adjust(wspace=0.1) 
fig.subplots_adjust(hspace=0.40) 
fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
fig.text(0.05,0.5, r"Helical fraction",ha='center',va='center',rotation='vertical') 
fig.subplots_adjust(left=0.1, bottom=0.1,right=0.80,top=0.9)

fig2, axarr2 = plt.subplots(figRows, figCols, sharex='none',sharey='row', figsize=(10.5,6.5)) 
fig2.subplots_adjust(wspace=0.1) 
fig2.subplots_adjust(hspace=0.40) 
fig2.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
fig2.text(0.05,0.5, r"Avg helical content",ha='center',va='center',rotation='vertical') 
fig2.subplots_adjust(left=0.1, bottom=0.1,right=0.80,top=0.9)

figD, axarrD = plt.subplots(figRows, figCols, sharex='none',sharey='row', figsize=(10.5,6.5)) 
figD.subplots_adjust(wspace=0.1) 
figD.subplots_adjust(hspace=0.40) 
figD.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
figD.text(0.05,0.5, r"Avg helical content",ha='center',va='center',rotation='vertical') 
figD.subplots_adjust(left=0.1, bottom=0.1,right=0.80,top=0.9)

blue = mpatches.Patch([],[],color='b',label=r"$\alpha$-helix")
green= mpatches.Patch([],[],color='g',label=r"$3_{10}$-helix")
red  = mpatches.Patch([],[],color='r',label=r"unfolded")

xmax=50
for row,solvent in enumerate(['sam','single_bound_sam','not_bound_sam']): 
    for col,state in enumerate(['folded','unfolded']) : 
        indFig = plt.figure() 
        #indAx = indFig.add_axes([0.1, 0.1, 0.7, 0.8])
        indFig.text(0.5,0.03, "Time (ns)", ha='center', va='center') 
        indFig.text(0.05,0.5, r"Helical fraction",ha='center',va='center',rotation='vertical') 
        indAx = indFig.add_subplot(111) 
        indFig.subplots_adjust(left=0.15, bottom=0.15,right=0.95,top=0.9)
        #indAx.set_xlabel("Time (ns)", ha='center', va='center') 
        #indAx.set_ylabel(r"Helical fraction",ha='center',va='center',rotation='vertical') 
        ax = axarr[row,col]
        ax2 = axarr2[row,col]
        axD = axarrD[row,col]
        datafile = '%s/%s/%s'%(solvent,state,inFile) 
        try : 
            data = np.genfromtxt(datafile) 
        except IOError : 
            print "No file found for %s/%s"%(solvent,state) 
            continue 

        if not endTime < 0 : 
            data = data[:endTime]
        frameNumber = len(data)
        assert frameNumber > 0 

        binNumber = frameNumber/binSize

        y1 = data[:,1] # Helical
        y2 = data[:,2] # B character
        y3 = data[:,3] # unfolded
        
        y1 /= 18 ##18 residues, so now each is fraction of total peptide
        y2 /= 18 
        y3 /= 18 

        ##Chop off data that does not fall into a bin
        while len(y1) %binSize !=0 : 
            y1 = y1[:-1] 
            y2 = y2[:-1] 
            y3 = y3[:-1] 
        assert len(y1) == len(y2) == len(y3) 

        y1binned = np.mean(y1.reshape(-1,binSize),axis=1) 
        y2binned = np.mean(y2.reshape(-1,binSize),axis=1) 
        y3binned = np.mean(y3.reshape(-1,binSize),axis=1) 

        y1average = np.cumsum(y1binned, dtype=float) 
        y2average = np.cumsum(y2binned, dtype=float) 
        y3average = np.cumsum(y3binned, dtype=float) 

        y1firstD = np.zeros_like(y1average) 
        y2firstD = np.zeros_like(y2average) 
        y3firstD = np.zeros_like(y3average) 

        for i in range(1,len(y1average)+1 )  : 
            y1average[i-1] /= i 
            y2average[i-1] /= i 
            y3average[i-1] /= i 

        for i in range(1,len(y1average) )  : 
            y1firstD[i] = y1average[i] - y1average[i-1]

        #rc_file('%s/rc_files/%s'%(projectDir,rcFile) ) 
        
        x = data[:,0] 
        x = np.linspace(0,np.max(data[:,0]),binNumber) 
        x /= 1000 
        if x.max() > xmax : 
            xmax = x.max()

        ax2.plot(x,y1average,color='b') 
        ax2.plot(x,y2average,color='g') 
        ax2.plot(x,y3average,color='r') 

        ax2.axhline(y1average[-1],color='b',linestyle='--',linewidth=1) 
        ax2.axhline(y2average[-1],color='g',linestyle='--',linewidth=1) 
        ax2.axhline(y3average[-1],color='r',linestyle='--',linewidth=1) 

        axD.plot(x,y1firstD,color='r') 
        axD.set_ylim(-0.0001,0.0001) 
        axD.axhline(0,color='k',linestyle='--',linewidth=.5) 

        ax2.scatter(x,y1binned,c='b',s=0.25,alpha=0.5)  
        ax2.scatter(x,y2binned,c='g',s=0.25,alpha=0.5)  
        ax2.scatter(x,y3binned,c='r',s=0.25,alpha=0.5)  

        ax.fill_between(x,0                ,y1binned                  ,facecolor='b',linewidth=0.0, edgecolor='none')
        ax.fill_between(x,y1binned         ,y1binned+y2binned         ,facecolor='g',linewidth=0.0, edgecolor='none')
        ax.fill_between(x,y1binned+y2binned,y1binned+y2binned+y3binned,facecolor='r',linewidth=0.0, edgecolor='none') 

        indAx.fill_between(x,0                ,y3binned         ,facecolor='r',linewidth=0.0, edgecolor='none')
        indAx.fill_between(x,y3binned         ,y3binned+y2binned,facecolor='g',linewidth=0.0, edgecolor='none') 
        indAx.fill_between(x,y3binned+y2binned,1                ,facecolor='b',linewidth=0.0, edgecolor='none') 
        
        ax.set_xlim(0,np.max(x)) 
        ax.set_ylim(0,1)
        ax.set_title("%s %s"%(state,solvent)) 

        indAx.set_xlim(0,xmax) 
        indAx.set_ylim(0,1)
        indAx.set_title("%s %s"%(state,solvent)) 
    
        if state == 'unfolded' : 
#            indFig.legend( [blue,green,red],[r"$\alpha$-helix",r"$3_{10}$-helix",r"unfolded"], 
#            loc = 'center', bbox_to_anchor=(0.90, 0.5),
#            fontsize='medium') 
            indAx.legend( [blue,green,red],[r"$\alpha$-helix",r"$3_{10}$-helix",r"unfolded"], 
            loc = 4,
            fontsize='medium') 
        indFig.savefig("%s/helicity_%s_%s.png"%(saveDir,state,solvent), format='png')
        plt.close(indFig) 
        print "%9s %9s plot complete"%(state,solvent)  

#fig.legend(handles=[blue, green, red],loc=4)
fig.legend( [blue,green,red],[r"helix",r"$\beta$",r"unfolded"], 
    loc = 'center', bbox_to_anchor=(0.90, 0.5),
    fontsize='medium') 

fig.savefig(outname, format='png')
fig2.savefig(outname2,format='png') 
figD.savefig(outnameD,format='png') 


plt.close() 
print "Combined plot complete." 




