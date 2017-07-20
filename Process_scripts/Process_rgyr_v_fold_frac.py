#!/usr/bin/env python

rcFile='paper.rc'
inFile ='rgyr/gyrate.xvg'
inFile2='dssp/helen.nrt'
saveDir='figures'
outname='%s/combined_rgyr_v_fold_frac.pdf'%saveDir

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm 
import glob
import os
import sys
from matplotlib import rc_file

figCols=2
figRows=3
xLabel=r"Radius of gyration (nm)"
yLabel=r"Helical fraction"

def Usage():
    print "Usage: %s <binSize>"%(sys.argv[0])

binSize=100
try : 
    binSize = int(sys.argv[1]) 
except : 
    Usage() 
    print "No bin size indicated, defaulting to %i"%binSize

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir) 

fig, axarr = plt.subplots(figRows, figCols, sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.1,hspace=0.25,left=0.1,bottom=0.1,right=0.80,top=0.9) 
fig.text(0.5,0.05, xLabel, ha='center', va='center') 
fig.text(0.05,0.5, yLabel, ha='center', va='center',rotation='vertical') 

#legend

for row, solvent in enumerate(['water','sam','tert']) : 
    for col,state in enumerate(['folded','unfolded']) : 
        indFig = plt.figure() 
        indAx = indFig.add_axes([0.1,0.1,0.7, 0.8]) 
        indFig.text(0.5,0.05, xLabel, ha='center', va='center') 
        indFig.text(0.05,0.5, yLabel, ha='center', va='center',rotation='vertical') 
        outnameInd = "%s/rgyr_v_fold_frac_%s_%s.pdf"%(saveDir,solvent,state) 

        ax = axarr[row,col]

        datafile = "%s_%s/%s"%(state,solvent,inFile)  
        datafile2= "%s_%s/%s"%(state,solvent,inFile2) 
        
        try : 
            data = np.genfromtxt(datafile,skip_header=25) 
            data2 = np.genfromtxt(datafile2) 
        except IOError : 
            print "No file found for %s %s"%(state,solvent) 
            break  
        except ValueError : 
            print "Tyring again with one less line..." 
            data = np.genfromtxt(datafile,skip_header=25,skip_footer=1) 
            data2 = np.genfromtxt(datafile2,skip_footer=1)  

        x = data[:-1,1] 
        y = (data2[:-1,1] + data2[:-1,2]) / 18

        while not len(x) == len(y) : 
            if len(x) > len(y) : 
                x = x[:-1]
            else : 
                y = y[:-1] 
        assert len(x) == len(y) 
        while len(x) % binSize !=0 : 
            x = x[:-1]
            y = y[:-1]
        assert len(x) % binSize == 0
        assert len(y) % binSize == 0
    
        ##Reshape into window averaged array
        xs = np.mean(x.reshape(-1,binSize),axis=1) 
        ys = np.mean(y.reshape(-1,binSize),axis=1) 
        
        colors = cm.brg(np.linspace(0,1,len(ys)) ) 
        for plot in ax, indAx : 
            plot.plot(xs,ys,color='k',alpha=0.5,zorder=1 ) 
            for x, y, c in zip(xs,ys, colors) : 
                plot.scatter(x ,y,color=c,edgecolor='none',s=40,alpha=1,zorder=2) 
            plot.set_xlim(0.8, 1.8) 
            plot.set_ylim(0,1.0) 
        #indAx.set_xlabel(xLabel) 
        #indAx.set_ylabel(yLabel,rotation='vertical') 

        indFig.savefig(outnameInd, format='pdf')
        plt.close(indFig)


fig.savefig(outname, format='pdf')
plt.close() 

