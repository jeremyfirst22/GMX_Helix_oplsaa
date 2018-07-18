#!/usr/bin/env python

rcFile='paper.rc'
inFile ='rgyr/gyrate.xvg'
inFile2='rmsd/rmsd.xvg'
saveDir='figures'
outname='%s/combined_rgyr_v_rmsd.pdf'%saveDir

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
yLabel=r"RMSD (nm)"

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
fig.subplots_adjust(wspace=0.1,hspace=0.25,left=0.1,bottom=0.1,right=0.8,top=0.9) 
fig.text(0.5,0.05, xLabel, ha='center', va='center') 
fig.text(0.05,0.5, yLabel, ha='center', va='center',rotation='vertical') 

#legend

for row, solvent in enumerate(['water','tert','sam']) : 
    for col,state in enumerate(['folded','unfolded']) : 
        indFig = plt.figure() 
        indAx = indFig.add_axes([0.1,0.1,0.9, 0.8]) 
        indFig.text(0.5,0.05, xLabel, ha='center', va='center') 
        indFig.text(0.05,0.5, yLabel, ha='center', va='center',rotation='vertical') 
        outnameInd = "%s/rgyr_v_rmsd_%s_%s.pdf"%(saveDir,solvent,state) 

        ax = axarr[row,col]

        datafile = "%s/%s/%s"%(solvent,state,inFile)  
        datafile2= "%s/%s/%s"%(solvent,state,inFile2) 
        
        try : 
            headlines = 0 
            with open(datafile) as f : 
                lines = f.readlines() 
                for line in lines : 
                    if line.startswith('#') or line.startswith('@') : 
                        headlines += 1 
                    else : 
                        break 
            data = np.genfromtxt(datafile,skip_header=headlines) 

            headlines = 0 
            with open(datafile2) as f : 
                lines = f.readlines() 
                for line in lines : 
                    if line.startswith('#') or line.startswith('@') : 
                        headlines += 1 
                    else : 
                        break 
            data2 = np.genfromtxt(datafile2,skip_header=headlines) 
        except IOError : 
            print "No file found for %s %s"%(state,solvent) 
            break  

        x = data[:-1,1] 
        y = data2[:-1,1]

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
        
        #colors = cm.brg(np.linspace(0,1,len(ys)) ) 
        cm = plt.cm.get_cmap('brg') 
        z = np.arange(0,len(ys) ) 
        z = z*binSize*4/1000
        for plot in ax, indAx : 
            plot.plot(xs,ys,color='k',alpha=0.5,zorder=1 ) 
            plot.set_xlim(0.8, 1.8) 
            plot.set_ylim(0,1.0) 
            plot.set_title('%s %s'%(state, solvent)) 
        sc = indAx.scatter(xs ,ys,c=z,edgecolor='none',s=40,alpha=1,zorder=2,vmin=0,vmax=max(z)+1) #frames -> ns
        ax.scatter(xs,ys,c=z,edgecolor='none',s=20, alpha=1, zorder=2, vmin=0, vmax=max(z)+1) 
        cbarInd = plt.colorbar(sc) 
        cbarInd.set_label('Time (ns)',rotation='vertical') 
        #indFig.colorbar(indAx) 
        #indAx.set_xlabel(xLabel) 
        #indAx.set_ylabel(yLabel,rotation='vertical') 

        indFig.savefig(outnameInd, format='pdf')
        plt.close(indFig)

cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
cbar = fig.colorbar(sc, cax=cbar_ax) 
cbar.set_label('Time (ns)',rotation='vertical') 
fig.savefig(outname, format='pdf')
plt.close() 

