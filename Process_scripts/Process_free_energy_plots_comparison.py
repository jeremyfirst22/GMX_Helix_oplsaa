#!/usr/bin/env python

rcFile='paper.rc'
inFile ='gibbs.dat'
saveDir='figures'
outname='%s/comparison_free_energy_plots.png'%(saveDir) 

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
figRows=2

rc_file("rc_files/%s"%rcFile) 

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir) 

colormap= plt.cm.get_cmap('jet')
figTA, axTA = plt.subplots(1,1)
Z = [[0,0],[0,0]]                                 # Apparently you can't set a colorbar without a mappable. This is hack
numLevels = 25          ##Number of levels to show in color bar, should be equal to number of levels calculated. 
vmin, vmax = -23.9, 0.0 ##min, max of dG in kJ/mol
levels = np.linspace(vmin,vmax,numLevels) #   around that. Taken from :
CS3 = axTA.contourf(Z, levels, cmap=colormap)     #   https://stackoverflow.com/questions/8342549/matplotlib-add-colorbar-to-a-sequence-of-line-plots
plt.close(figTA)


left, right = 0.150, 0.725
bottom, top =  0.10,0.98
hspace, wspace = 0.15,0.15

fig, axarr = plt.subplots(figRows, figCols, sharex='all',sharey='all', figsize=(3.25,3.0))
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top)
fig.subplots_adjust(wspace=wspace,hspace=hspace)

fig.text((right-left)/2+left,0,           r"Radius of gyration, $R_g$ ($\AA$)", ha='center', va='bottom')
fig.text(0.01,(top-bottom)/2+bottom,         r"RMSD ($\AA$)",ha='left',va='center',rotation='vertical')

fig.text(right,top-(top-bottom-hspace/2)/4 ,          r"Ordered SAM",ha='left',va='center',rotation=270)
fig.text(right,(top-bottom-hspace/2)/4+bottom, r"Disordered SAM",ha='left',va='center',rotation=270)

indexToLetter = {
        1:'A',
        2:'B',
        3:'C',
        }

index = 0 
for row, solvent in enumerate(['ordered', 'disordered']) : 
    xMin, yMin = -1, -1 
    xMax, yMax = 1e10, 1e10
    index += 1 

    ax = axarr[row]

    datafile = "sham_comparison/%s/%s"%(solvent,inFile)  

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
    im = ax.pcolor(x,y, data, cmap=colormap,vmin=vmin, vmax=vmax) 

    minY, minX = np.unravel_index(np.nanargmin(data), np.shape(data))
    ax.axhline(y[minY],linestyle='--', color='k',linewidth=0.5)
    ax.axvline(x[minX],linestyle='--', color='k',linewidth=0.5)

    ax.text(0.97,0.95,r"\textsf{%c}"%indexToLetter[index],va='top',ha='right',transform=ax.transAxes,fontsize=12)

cbarHeight = (top-bottom)*0.75 
cbarY = bottom + (top-bottom-cbarHeight)/2

cbar_ax = fig.add_axes([0.85, cbarY, 0.03, cbarHeight])
cbar = fig.colorbar(CS3, cax=cbar_ax,ticks=np.arange(vmax,vmin,-2)) 

fig.text(0.85, cbarY + cbarHeight/2 ,r"$\Delta G$ (kJ mol$^{-1}$)",ha='right',va='center',rotation='vertical') 
ax.set_xlim([7.5,17.5])
ax.set_ylim([1.50, 11]) 

fig.savefig(outname, format='png')
plt.close() 

