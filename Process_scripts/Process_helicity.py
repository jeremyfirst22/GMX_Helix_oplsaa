#!/usr/bin/env python

projectDir='/Users/jeremyfirst/GMX_Helix_oplsaa'
#rcFile='poster.rc'
rcFile='paper.rc'
inFile='dssp/helen.nrt'
saveDir='figures'
outname='%s/combined_helicity.png'%saveDir
outname2='%s/combined_average_helicity.png'%saveDir
outnameD='%s/combined_firstD_helicity.png'%saveDir

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

left, right = 0.08, 0.8 
top, bottom = 0.95, 0.1
hspace, wspace = 0.15,0.1

fig, axarr = plt.subplots(figRows, figCols, sharex='col',sharey='row', figsize=(6.5,3.5)) 
fig.subplots_adjust(wspace=wspace) 
fig.subplots_adjust(hspace=hspace) 
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top)

fig.text((right-left)/2+left,0.03,           r"Time (ns)", ha='center', va='center') 
fig.text(0.02,(top-bottom)/2+bottom,         r"Helical fraction",ha='center',va='center',rotation='vertical') 
fig.text((right-left)/4+left,top,            r"Starting with folded",ha='center',va='bottom') 
fig.text(right-(right-left)/4,top,           r"Starting with unfolded",ha='center',va='bottom') 
fig.text(right,top-(top-bottom)/6 ,          r"H$_2$O",ha='left',va='center',rotation=270) 
fig.text(right,top-(top-bottom)/2,           r"2:1 H$_2$O:$t$-BuOH",ha='left',va='center',rotation=270) 
fig.text(right,(top-bottom-hspace)/6+bottom, r"SAM surface",ha='left',va='center',rotation=270) 

fig2, axarr2 = plt.subplots(figRows, figCols, sharex='col',sharey='row', figsize=(6.5,3.5)) 
fig2.subplots_adjust(wspace=wspace) 
fig2.subplots_adjust(hspace=hspace) 
fig2.subplots_adjust(left=left, bottom=bottom,right=right,top=top)

fig2.text((right-left)/2+left,0.03,           r"Time (ns)", ha='center', va='center') 
fig2.text(0.02,(top-bottom)/2+bottom,         r"Average helical fraction",ha='center',va='center',rotation='vertical') 
fig2.text((right-left)/4+left,top,            r"Starting with folded",ha='center',va='bottom') 
fig2.text(right-(right-left)/4,top,           r"Starting with unfolded",ha='center',va='bottom') 
fig2.text(right,top-(top-bottom)/6 ,          r"H$_2$O",ha='left',va='center',rotation=270) 
fig2.text(right,top-(top-bottom)/2,           r"2:1 H$_2$O:$t$-BuOH",ha='left',va='center',rotation=270) 
fig2.text(right,(top-bottom-hspace)/6+bottom, r"SAM surface",ha='left',va='center',rotation=270) 


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
for row,solvent in enumerate(['water','tert','not_bound_sam']): 
    y1Total,y2Total,y3Total = 0,0,0
    for col,state in enumerate(['folded','unfolded']) : 
        continue 
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

        if solvent == 'water' : 
            equilTime = 600 ##ns
        elif solvent == 'tert' : 
            if state == 'folded' : 
                equilTime = 150 
            if state == 'unfolded' : 
                equilTime = 250 
        elif solvent == 'not_bound_sam' : 
            equilTime = 50 

        equilTime *= 1000 ##ns -> ps 
        equilTime /= 4 ##ps -> frames 

        #binNumber = (frameNumber-equilTime)/binSize
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

        y1average = np.cumsum(y1binned[equilTime/binSize:], dtype=float) 
        y2average = np.cumsum(y2binned[equilTime/binSize:], dtype=float) 
        y3average = np.cumsum(y3binned[equilTime/binSize:], dtype=float) 

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
        
        #x = data[equilTime:,0] 
        #x = np.linspace((equilTime*4),np.max(data[:,0]),binNumber) 
        x = np.linspace(0,np.max(data[:,0]),binNumber) 
        x /= 1000 
        if x.max() > xmax : 
            xmax = x.max()

        ax2.plot(x[equilTime/binSize:],y1average,color='b') 
        ax2.plot(x[equilTime/binSize:],y2average,color='g') 
        ax2.plot(x[equilTime/binSize:],y3average,color='y') 
        
        ax2.axvline(equilTime*4/1000,color='k') 

        ax2.axhline(y1average[-1],color='b',linestyle='--',linewidth=1) 
        ax2.axhline(y2average[-1],color='g',linestyle='--',linewidth=1) 
        ax2.axhline(y3average[-1],color='y',linestyle='--',linewidth=1) 

        #axD.plot(x,y1firstD,color='r') 
        #axD.set_ylim(-0.0001,0.0001) 
        #axD.axhline(0,color='k',linestyle='--',linewidth=.5) 

        ax2.scatter(x,y1binned,c='b',s=0.15,alpha=0.5)  
        ax2.scatter(x,y2binned,c='g',s=0.15,alpha=0.5)  
        ax2.scatter(x,y3binned,c='r',s=0.15,alpha=0.5)  

        ax.fill_between(x,0                ,y1binned                  ,facecolor='b',linewidth=0.0, edgecolor='none')
        ax.fill_between(x,y1binned         ,y1binned+y2binned         ,facecolor='g',linewidth=0.0, edgecolor='none')
        ax.fill_between(x,y1binned+y2binned,y1binned+y2binned+y3binned,facecolor='r',linewidth=0.0, edgecolor='none') 

        ax.set_xlim(0,np.max(x)) 
        ax.set_ylim(0,1)
        
        #ax.set_title('%s %s'%(state,solvent.replace('_',' ')) )

        #print "%9s %9s plot complete"%(state,solvent)  

        y1Total += y1average[-1]
        y2Total += y2average[-1]
        y3Total += y3average[-1]

        print "%10s\t%10s\tHelix: %0.3f\tBeta: %0.3f\tUnfolded: %0.3f"%(state,solvent,y1average[-1],y2average[-1],y3average[-1]) 
    print "%10s\t%10s\tHelix: %0.3f\tBeta: %0.3f\tUnfolded: %0.3f"%("Total: ",solvent,y1Total/2, y2Total/2, y3Total/2)

#fig.legend(handles=[blue, green, red],loc=4)
fig.legend( [blue,green,red],[r"helix",r"$\beta$",r"unfolded"], 
    loc=(0.85,0.47),
    fontsize='medium') 

fig.savefig(outname, format='png')
fig2.savefig(outname2,format='png') 
figD.savefig(outnameD,format='png') 


plt.close() 
print "Combined plot complete." 




