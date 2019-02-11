#!/usr/bin/env python

projectDir='/Users/jeremyfirst/GMX_Helix_oplsaa'
#rcFile='poster.rc'
rcFile='paper.rc'
dataFile='dssp/scount.xvg'
saveDir='figures'
outname='%s/bounds_helicity.png'%saveDir
outname2='%s/bounds_average_helicity.png'%saveDir

figCols=2
figRows=2

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import glob
import os
import sys
from matplotlib import rc_file
import matplotlib as mpl

rc_file('rc_files/%s'%rcFile) 

#mpl.rcParams['figure.figsize'] = 5,3.8

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
top, bottom = 0.93, 0.13
hspace, wspace = 0.15,0.1

fig, axarr = plt.subplots(figRows, figCols, sharex='col',sharey='row', figsize=(6.5,2.3)) 
fig.subplots_adjust(wspace=wspace) 
fig.subplots_adjust(hspace=hspace) 
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top)

fig.text((right-left)/2+left,0.03,           r"Time (ns)", ha='center', va='center') 
fig.text(0.02,(top-bottom)/2+bottom,         r"Conformational fraction",ha='center',va='center',rotation='vertical') 
fig.text((right-left)/4+left,top,            r"Starting with folded",ha='center',va='bottom') 
fig.text(right-(right-left)/4,top,           r"Starting with unfolded",ha='center',va='bottom') 
fig.text(right,(top - bottom)/4+bottom,      r"1 restraint",ha='left',va='center',rotation=270) 
fig.text(right,top-(top-bottom)/4,           r"2 restraints",ha='left',va='center',rotation=270) 

fig2, axarr2 = plt.subplots(figRows, figCols, sharex='col',sharey='row', figsize=(6.5,2.3)) 
fig2.subplots_adjust(wspace=wspace) 
fig2.subplots_adjust(hspace=hspace) 
fig2.subplots_adjust(left=left, bottom=bottom,right=right,top=top)

fig2.text((right-left)/2+left,0.03,           r"Time (ns)", ha='center', va='center') 
fig2.text(0.02,(top-bottom)/2+bottom,         r"Average helical fraction",ha='center',va='center',rotation='vertical') 
fig2.text((right-left)/4+left,top,            r"Starting with folded",ha='center',va='bottom') 
fig2.text(right-(right-left)/4,top,           r"Starting with unfolded",ha='center',va='bottom') 
fig2.text(right,(top - bottom)/4+bottom,      r"1 restraint ",ha='left',va='center',rotation=270) 
fig2.text(right,top-(top-bottom)/4,           r"2 restraints",ha='left',va='center',rotation=270) 

blue = mpatches.Patch([],[],color='b',label=r"Helical")
royalblue = mpatches.Patch([],[],color='royalblue',label=r"Turn")
green= mpatches.Patch([],[],color='g',label=r"Strand") 
red  = mpatches.Patch([],[],color='r',label=r"Unfolded")

xmax=50
for row,solvent in enumerate(['sam','single_bound_sam']) : 
    y1Total,y2Total,y3Total,y1dTotal = 0,0,0,0
    for col,state in enumerate(['folded','unfolded']) : 
        try : 
            del alpha,three,five,turn,bsheet,bbridge,bend,coil
        except : 
            pass

        inFile="%s/%s/%s"%(solvent,state,dataFile) 

        ax = axarr[row,col]
        ax2 = axarr2[row,col]

        headlines=0 
        with open(inFile) as f : 
            for line in f.readlines() : 
                if line.startswith('#') or line.startswith('@') : 
                    headlines += 1 
                    if 'A-Helix' in line :
                        alpha=line.split()[1][1:]
                    if '3-Helix' in line :
                        three=line.split()[1][1:]
                    if '5-Helix' in line :
                        five =line.split()[1][1:]
                    if 'Turn' in line :
                        turn =line.split()[1][1:]
                    if 'B-Sheet' in line :
                        bsheet = line.split()[1][1:]
                    if 'B-Bridge' in line :
                        bbridge= line.split()[1][1:]
                    if 'Bend' in line :
                        bend = line.split()[1][1:]
                    if 'Coil' in line :
                        coil = line.split()[1][1:]
                else : 
                    break 

        data = np.genfromtxt(inFile,skip_header=headlines,skip_footer=2) 

        try :
            y1a = data[:,int(alpha)+1]
        except NameError :
            y1a = np.zeros(len(data[:,0]))
        try :
            y1b = data[:,int(three)+1]
        except NameError :
            y1b = np.zeros(len(data[:,0]))
        try :
            y1c = data[:,int(five)+1]
        except NameError :
            y1c = np.zeros(len(data[:,0]))
        try :
            y1d = data[:,int(turn)+1]
        except NameError :
            y1d = np.zeros(len(data[:,0]))
        y1 = y1a + y1b + y1c #+ y1d
        
        try :
            y2a = data[:,int(bsheet)+1]
        except NameError :
            y2a = np.zeros(len(data[:,0]))
        try :
            y2b = data[:,int(bbridge)+1]
        except NameError :
            y2b = np.zeros(len(data[:,0]))
        try :
            y2c = data[:,int(bend)+1]
        except NameError :
            y2c = np.zeros(len(data[:,0]))
        y2 = y2a + y2b + y2c 
        
        try : 
            y3 = data[:,int(coil)+1]
        except NameError : 
            y3 = np.zeros(len(data[:,0])) 
        y3 = 20 - y1 - y2 - y1d
        
        if not endTime < 0 : 
            data = data[:endTime]
        frameNumber = len(data)
        assert frameNumber > 0 

        if solvent == 'water' : 
            equilTime = 600 ##ns
        else : 
            equilTime = 150 ##ns 

        equilTime *= 1000 ##ns -> ps 
        equilTime /= 4 ##ps -> frames 

        binNumber = frameNumber/binSize
        
        y1 /= 20 ##20 residues, so now each is fraction of total peptide
        y2 /= 20 
        y3 /= 20 
        y1d/= 20 

        #Chop off data that does not fall into a bin
        while len(y1) %binSize !=0 : 
            y1 = y1[:-1] 
            y1d= y1d[:-1] 
            y2 = y2[:-1] 
            y3 = y3[:-1] 
        assert len(y1) == len(y2) == len(y3) 

        y1binned = np.mean(y1.reshape(-1,binSize),axis=1) 
        y2binned = np.mean(y2.reshape(-1,binSize),axis=1) 
        y3binned = np.mean(y3.reshape(-1,binSize),axis=1) 
        y1dbinned = np.mean(y1d.reshape(-1,binSize),axis=1) 

        y1average = np.cumsum(y1binned[equilTime/binSize:], dtype=float) 
        y2average = np.cumsum(y2binned[equilTime/binSize:], dtype=float) 
        y3average = np.cumsum(y3binned[equilTime/binSize:], dtype=float) 
        y1daverage = np.cumsum(y1dbinned[equilTime/binSize:], dtype=float) 

        for i in range(1,len(y1average)+1 )  : 
            y1average[i-1] /= i 
            y2average[i-1] /= i 
            y3average[i-1] /= i 
            y1daverage[i-1] /= i 

        x = np.linspace(0,np.max(data[:,0]),binNumber) 
        x /= 1000 
        if x.max() > xmax : 
            xmax = x.max()

        ax2.plot(x[equilTime/binSize:],y1average,color='b') 
        ax2.plot(x[equilTime/binSize:],y2average,color='g') 
        ax2.plot(x[equilTime/binSize:],y3average,color='y') 
        ax2.plot(x[equilTime/binSize:],y1daverage,color='b') 
        
        ax2.axvline(equilTime*4/1000,color='k',linestyle='--') 
        ax.axvline(equilTime*4/1000,color='k',linestyle='--') 

        ax2.axhline(y1average[-1],color='b',linestyle='--',linewidth=1) 
        ax2.axhline(y2average[-1],color='g',linestyle='--',linewidth=1) 
        ax2.axhline(y3average[-1],color='y',linestyle='--',linewidth=1) 
        ax2.axhline(y1daverage[-1],color='b',linestyle='--',linewidth=1) 

        ax2.scatter(x,y1binned,c='b',s=0.15,alpha=0.5)  
        ax2.scatter(x,y2binned,c='g',s=0.15,alpha=0.5)  
        ax2.scatter(x,y3binned,c='r',s=0.15,alpha=0.5)  
        ax2.scatter(x,y1dbinned,c='b',s=0.15,alpha=0.5)  

        ax.fill_between(x,0                ,y1binned                  ,facecolor='b',linewidth=0.0, edgecolor='none')
        ax.fill_between(x,y1binned         ,y1binned+y1dbinned        ,facecolor='royalblue',linewidth=0.0, edgecolor='none')
        ax.fill_between(x,y1binned+y1dbinned         ,y1binned+y1dbinned+y2binned         ,facecolor='g',linewidth=0.0, edgecolor='none')
        ax.fill_between(x,y1binned+y1dbinned+y2binned,y1binned+y1dbinned+y2binned+y3binned,facecolor='r',linewidth=0.0, edgecolor='none')

        ax.set_xlim(0,np.max(x)) 
        ax.set_ylim(0,1)
        #ax.set_title("%s %s"%(state,solvent.replace('_',' '))) 

        y1Total += y1average[-1]
        y2Total += y2average[-1]
        y3Total += y3average[-1]
        y1dTotal += y1daverage[-1]

        print "%10s\t%10s\tHelix: %0.3f\tTurn: %0.3f\tBeta: %0.3f\tUnfolded: %0.3f"%(state,solvent,y1average[-1],y1daverage[-1],y2average[-1],y3average[-1])
    print "%10s\t%10s\tHelix: %0.3f\tTurn: %0.3f\tBeta: %0.3f\tUnfolded: %0.3f"%("Total: ",solvent,y1Total/2,y1dTotal/2,y2Total/2, y3Total/2)

fig.legend( [blue,royalblue,green,red],[r"Helix","Turn",r"Strand",r"Unfolded"],
    loc = 'center', bbox_to_anchor=(0.90, 0.525),
    fontsize='medium')

fig.savefig(outname, format='png')
fig2.savefig(outname2,format='png') 

plt.close() 
print "Combined plot complete." 




