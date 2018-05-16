#!/usr/bin/env python

inFile='density/charge.xvg'
saveDir='figures'
outname='%s/charge.pdf'%saveDir

figCols=2
figRows=7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import glob
import os
import sys
from matplotlib import rc_file
import matplotlib as mpl
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

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

fig, axarr = plt.subplots(figRows, figCols, sharex='col',sharey='row',figsize=( 7,11)) 
fig.subplots_adjust(wspace=0.1) 
fig.subplots_adjust(hspace=0.25) 
fig.text(0.5,0.05, r"Distance ($\rm{\AA}$)", ha='center', va='center') 
fig.text(0.05,0.5, r"Denisty (kg/m$^3$)",ha='center',va='center',rotation='vertical') 
fig.subplots_adjust(left=0.15, bottom=0.1,right=0.80,top=0.9)

for row, salt in enumerate(np.arange(0,175,25)) :
    for col,state in enumerate(['folded','unfolded']) : 
        ax = axarr[row,col]
        datafile = 'sam_%s/%s/%s'%(salt,state,inFile) 
        try : 
            headerlines = 0 
            with open(datafile) as f: 
                fileLines = f.readlines() 
                for line in fileLines : 
                    if line.startswith('#') or line.startswith('@') :
                        headerlines+=1
                        continue 
                    else : 
                        break        
            data = np.genfromtxt(datafile,skip_header=headerlines) 
        except IOError : 
            print "No file found for %s/%s"%(salt,state) 
            continue 

        ax.plot(data[:,0], data[:,1],linewidth=0.5) 
        ax.plot(data[:,0], data[:,2],linewidth=0.5) 
        ax.plot(data[:,0], data[:,3],linewidth=0.5) 
        ax.plot(data[:,0], data[:,4],linewidth=0.5) 
        if salt != 0 : 
            ax.plot(data[:,0], data[:,5],linewidth=0.5) 

        inset = inset_axes(ax,
                    width="50%", # width = 30% of parent_bbox
                    height="55%", # height : 1 inch
                    loc='center right') 
        inset.plot(data[:,0], data[:,1],linewidth=0.5) 
        inset.plot(data[:,0], data[:,2],linewidth=0.5) 
        inset.plot(data[:,0], data[:,3],linewidth=0.5) 
        inset.plot(data[:,0], data[:,4],linewidth=0.5) 
        if salt != 0 : 
            inset.plot(data[:,0], data[:,5],linewidth=0.5) 

#        ax.set_xlim(0,8) 
#        ax.set_ylim(0,1200) 

#        inset.set_xlim(2,8) 
#        inset.set_ylim(0,15) 


fig.savefig(outname, format='pdf')
plt.close()

print "Combined plot complete." 




