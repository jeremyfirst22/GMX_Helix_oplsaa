#!/usr/bin/env python

projectDir='/Users/jeremyfirst/GMX_Helix_oplsaa'
#rcFile='poster.rc'
rcFile='paper.rc'
inFile='rgyr/gyrate.xvg'
saveDir='figures'
#outname='%s/combined_helicity.pdf'%saveDir

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import glob
import os
import sys
from matplotlib import rc_file

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir) 

figRows, figCols = 1,3 

for row,solvent in enumerate(['water','tert','sam']): 
    #fig, axarr = plt.subplots(1,3,sharey='row',gridspec_kw = {'width_ratios':[3,3,1]}) 
    fig, axarr = plt.subplots(figRows, figCols, sharex='none',sharey='row', figsize=(10.5,6.5), gridspec_kw = {'width_ratios':[3,3,1]})
    fig.subplots_adjust(wspace=0.1)
    fig.subplots_adjust(hspace=0.40)
    fig.text(0.5,0.04, "Time (ns)", ha='center', va='center')
    fig.text(0.05,0.5, r"Helical fraction",ha='center',va='center',rotation='vertical')
    fig.subplots_adjust(left=0.1, bottom=0.1,right=0.80,top=0.9)
    fig.subplots_adjust(wspace=0) 
    for col,state in enumerate(['folded','unfolded']) : 
        ax = axarr[col]

        color = 'b'
        if col == 0 : 
            color = 'r'

        datafile = '%s/%s/%s'%(solvent,state,inFile) 
        try : 
            headlines = 0
            with open(datafile) as f: 
                for line in f.readlines() : 
                    if line.startswith('#') or line.startswith('@') : 
                        headlines += 1 
                    else : 
                        break 
            data = np.genfromtxt(datafile,skip_header=headlines) 
        except IOError : 
            print "No file found for %s/%s"%(solvent,state) 
            continue 
        data[:,0] /= 1000  ##ps -> ns 

        #data = data[25000:,:]
        
        #rc_file('%s/rc_files/%s'%(projectDir,rcFile) ) 
        
        ax.plot(data[:,0],data[:,1],color=color) 


        ax3 = axarr[2] 
        x,y = np.histogram(data[:,1], np.linspace(0.6,1.8,100),density=True) 
        y = y[:-1]

        ax3.plot(x,y,color=color)

        #ax.set_xlim(0,xmax) 
        ax.set_ylim(0.6,1.8) 
        #ax.set_title("%s %s"%(state,solvent)) 
    
    fig.savefig("%s/rgyr_%s_%s.png"%(saveDir,state,solvent), format='png')
    print "%9s %9s plot complete"%(state,solvent)  





