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

left, right = 0.25,0.95
bottom, top = 0.20,0.95

fig, ax = plt.subplots(1,1,figsize=(2.5,1.5))
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top)
fig.text((right-left)/2+left,0.00,           r"Height ($\AA$)", ha='center', va='bottom')
fig.text(0.01,(top-bottom)/2+bottom,         r"Density (kg m$^{-3}$)",ha='left',va='center',rotation='vertical')

for solvent in ['disordered_not_bound_sam','not_bound_sam'] :
    datafile = "%s/folded/order/density.xvg"%solvent

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

    data[:,0] *= 10 

    if solvent == 'not_bound_sam' : 
        label="Ordered SAM"
        color='k'
    else : 
        label="Disordered SAM" 
        color='r'

    ax.plot(data[:,0],data[:,1],label=label,color=color) 

ax.legend(loc=1)
ax.set_xlim([5,25]) 

fig.savefig('figures/density.png',format='png') 

