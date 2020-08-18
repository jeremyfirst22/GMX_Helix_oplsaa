#!/usr/bin/env python

rcFile='paper.rc'
inFile ='sham/gibbs.dat'
saveDir='figures'
outname='%s/free_energy_plots.png'%(saveDir) 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm 
from matplotlib import rc_file
import glob
import os
import sys
from matplotlib import rc_file
from adjustText import adjust_text

figCols=3
figRows=1

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

left, right = 0.07, 0.90
bottom, top =  0.20,0.90
hspace, wspace = 0.00,0.15

fig, axarr = plt.subplots(figRows, figCols, sharex='all',sharey='all', figsize=(6.50,1.5))
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top)
fig.subplots_adjust(wspace=wspace,hspace=hspace)

fig.text((right-left)/2+left,0,           r"Radius of gyration, $R_g$ ($\AA$)", ha='center', va='bottom')
fig.text(0.01,(top-bottom)/2+bottom,         r"RMSD ($\AA$)",ha='left',va='center',rotation='vertical')

fig.text((right-left-wspace/2)/6 +left          ,0.98, r"H$_2$O",ha='center',va='top') 
fig.text((right-left)/2          +left          ,0.98, r"2:1 $t$-BuOH/H$_2$O",ha='center',va='top') 
fig.text(right-(right-left-wspace/2)/6          ,0.98, r"SAM surface",ha='center',va='top') 

indexToLetter = {
        1:'A',
        2:'B',
        3:'C',
        }

index = 0 
for row, solvent in enumerate(['water','tert','not_bound_sam']) : 
    labels = [] 
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

    data -= np.max(data)    ##Set max dG (not sampled region) to zero. dG are relative to that. 

    data[data == 0] = 'nan' ##Replace dG = 0 with 'nan' so it doesn't show up on heat map. 

    x = np.linspace(xMin,xMax,np.shape(data)[0]) 
    x *= 10 ## nm -> A
    y = np.linspace(yMin,yMax,np.shape(data)[1]) 
    y *= 10 ## nm -> A
    
    print np.nanmin(data) 
    im = ax.pcolor(x,y, data, cmap=colormap,vmin=vmin, vmax=vmax) 

    minY, minX = np.unravel_index(np.nanargmin(data), np.shape(data))
    ax.axhline(y[minY],linestyle='--', color='k',linewidth=0.5) 
    ax.axvline(x[minX],linestyle='--', color='k',linewidth=0.5) 

    ##Label rgyr and rmsd of top clusters
    numStructures = 5  ##number of top clusters to label
    rmsdFile = "%s/rgyr_rmsd_cluster/rmsd.xvg"%(solvent) 
    try : 
        headlines = 0 
        with open(rmsdFile) as f : 
            for line in f.readlines() : 
                if line.startswith('#') or line.startswith('@') : 
                    headlines += 1 
                else : 
                    break 
        rmsdDat = np.genfromtxt(rmsdFile, skip_header=headlines)
    except IOError : 
        print "%s file not found. Skipping labeling."%rmsdFile 
        continue 
    rgyrFile = "%s/rgyr_rmsd_cluster/gyrate.xvg"%(solvent) 
    try : 
        headlines = 0 
        with open(rgyrFile) as f : 
            for line in f.readlines() : 
                if line.startswith('#') or line.startswith('@') : 
                    headlines += 1 
                else : 
                    break 
        rgyrDat = np.genfromtxt(rgyrFile, skip_header=headlines)
    except IOError : 
        print "%s file not found. Skipping labeling."%rgyrFile 
        continue 
    for i in range(numStructures) : 
        labels.append(ax.text(rgyrDat[i][1]*10,rmsdDat[i][1]*10, "%i"%(i+1), va='center', ha='center') )
        continue

    adjust_text(labels,ax=ax) #,arrowprops=dict(arrowstyle='-')  ) 

    for i in range(numStructures) : 
        ax.scatter(rgyrDat[i][1]*10,rmsdDat[i][1]*10, s=0.25,c='k')
        continue


    ax.text(0.97,0.95,r"\textsf{%c}"%indexToLetter[index],va='top',ha='right',transform=ax.transAxes,fontsize=12) 

leftbar = (1-right)/4+right+0.01
widthbar = (1-right)/6 

cbar_ax = fig.add_axes([leftbar, bottom, widthbar, (top-bottom)])
cbar = fig.colorbar(CS3, cax=cbar_ax,orientation='vertical') 
cbar.set_ticks(np.arange(-20,1,5)) 

fig.text(leftbar, (top-bottom)/2+bottom,r"$\Delta G$ (kJ mol$^{-1}$)",ha='right',va='center',rotation='vertical') 
ax.set_xlim([7.5,17.5])
ax.set_ylim([1.50, 11]) 

fig.savefig(outname, format='png')
plt.close() 

