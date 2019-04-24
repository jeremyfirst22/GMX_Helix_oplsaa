#!/usr/bin/env python

rcFile='paper.rc'
inFile ='sham/gibbs.dat'
saveDir='figures'
outname='%s/disorder_free_energy_plots.png'%(saveDir) 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm 
from matplotlib import rc_file
import glob
import os
import sys
from matplotlib import rc_file

figCols=1
figRows=3

rc_file("rc_files/%s"%rcFile) 

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir) 

left, right = 0.15, 0.65
bottom, top =  0.1,0.98
hspace, wspace = 0.15,0.15

fig, axarr = plt.subplots(figRows, figCols, sharex='all',sharey='all', figsize=(3.25,3.8))
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top)
fig.subplots_adjust(wspace=wspace,hspace=hspace)

fig.text((right-left)/2+left,0,           r"Radius of gyration R$_g$ ($\AA$)", ha='center', va='bottom')
fig.text(0.01,(top-bottom)/2+bottom,         r"RMSD ($\AA$)",ha='left',va='center',rotation='vertical')

fig.text(right,top-(top-bottom)/6 ,          r"2 restraints",ha='left',va='center',rotation=270)
fig.text(right,top-(top-bottom)/2,           r"1 restraint",ha='left',va='center',rotation=270)
fig.text(right,(top-bottom-hspace)/6+bottom, r"no restraints",ha='left',va='center',rotation=270)

indexToLetter = {
        1:'A',
        2:'B',
        3:'C',
        }

index = 0 
for row, solvent in enumerate(['disordered_sam','disordered_single_bound_sam','disordered_not_bound_sam']) : 
    xMin, yMin = -1, -1 
    xMax, yMax = 1e10, 1e10
    index += 1 

    ax = axarr[row]

    datafile = "%s/%s"%(solvent,inFile)  

    with open(datafile) as f : 
        lines = f.readlines() 
        for line in lines[:2] : 
            if 'x-range:' in line : 
                xMin = float(line.split()[2])
                xMax = float(line.split()[3])
            if 'y-range:' in line : 
                yMin = float(line.split()[2])
                yMax = float(line.split()[3])
    
    data = np.genfromtxt(datafile,skip_header=2) 

    data -= np.max(data) 

    data[data == 0] = 'nan'

    x = np.linspace(xMin,xMax,np.shape(data)[0]) 
    x *= 10 ## nm -> A
    y = np.linspace(yMin,yMax,np.shape(data)[1]) 
    y *= 10 ## nm -> A
    
    print np.nanmin(data) 
    im = ax.pcolor(x,y, data, cmap='Spectral_r',vmin=-22.5, vmax=0.0) 

cbar_ax = fig.add_axes([0.75, 0.15, 0.03, 0.75])
cbar = fig.colorbar(im, cax=cbar_ax) 
fig.text(0.75, 0.13         ,     r"$\Delta G$ (kJ/mol)",ha='left',va='top') 
ax.set_xlim([7.5,17.5])
ax.set_ylim([1.50, 11]) 

fig.savefig(outname, format='png')
plt.close() 

