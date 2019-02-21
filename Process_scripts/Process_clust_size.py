#!/usr/bin/env python 

projectDir='Users/jeremyfirst/GMX_Helix_oplsaa'
#rcFile='paper.rc'
rcFile='paper.rc'
saveDir='figures'

import numpy as np  
import matplotlib.pyplot as plt 
from matplotlib import rc_file
import glob as glob 
import os  
import sys 

rc_file('rc_files/%s'%rcFile) 
colormap= plt.cm.get_cmap('viridis')

top=os.getcwd()

colorDict = {'water':'k','tert':'b','not_bound_sam':'g'}

left, right = 0.13,0.95
bottom, top = 0.13,0.95

fig, ax = plt.subplots(1,1,figsize=(3.25,2))
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top)
fig.text((right-left)/2+left,0.00,           r"\# of clusters", ha='center', va='bottom')
fig.text(0.01,(top-bottom)/2+bottom,         r"\% of total ensemble",ha='left',va='center',rotation='vertical')

for solvent in ['tert','water','not_bound_sam' ] : 
    datafile = "%s/cluster/clust-size.xvg"%solvent

    if not os.path.isfile(datafile) : 
        print "%s file not found. Skipping."%solvent 
        continue 
    
    headlines = 0 
    with open(datafile) as f : 
        for line in f.readlines() : 
            if line.startswith('#') or line.startswith('@') : 
                headlines += 1 
            else : 
                break 
    data = np.genfromtxt(datafile,skip_header=headlines) 

    if solvent == 'tert' : 
        label = r"2:1 \emph{t}-BuOH:H$_2$O"
    elif solvent == 'water' : 
        label = r"H$_2$O"
    elif solvent == 'not_bound_sam' : 
        label = r"SAM surface"
    ax.plot(data[:,0],np.cumsum(data[:,1]/np.sum(data[:,1]) * 100) ,color=colorDict[solvent],label=label) 

#ax.set_yscale('log') 
ax.set_xscale('log') 
ax.set_xlim([1,1000]) 
ax.set_ylim([0,100]) 
ax.legend(loc=4)

fig.savefig('figures/clust-size.png',format='png') 

