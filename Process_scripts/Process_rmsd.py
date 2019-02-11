#!/usr/bin/env python

projectDir='/Users/jeremyfirst/GMX_Helix_oplsaa'
#rcFile='poster.rc'
rcFile='paper.rc'
inFile='rmsd/rmsd.xvg'
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

fig = plt.figure() 
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
fig.text(0.5,0.05, "Time (ns)", ha='center', va='center') 
fig.text(0.05,0.5, "RMSD (A )",ha='center',va='center',rotation='vertical') 

fig2= plt.figure() 
ax2= fig2.add_axes([0.1, 0.1, 0.8, 0.8])
fig2.text(0.5,0.05, "Time (ns)", ha='center', va='center') 
fig2.text(0.05,0.5, "RMSD (A )",ha='center',va='center',rotation='vertical') 

for row,solvent in enumerate(['water','tert','not_bound_sam']): 
    if solvent == 'water' :
        equilTime = 600 ##ns
    else :
        equilTime = 150 ##ns
    equilTime *= 1000 # ns -> ps
    equilTime /= 4 ##ps -> frames
    for col,state in enumerate(['folded','unfolded']) : 
        if solvent == 'not_bound_sam' : 
            color = 'g'
        elif solvent == 'tert' : 
            color = 'b' 
        else : 
            color = 'k' 

        if state == 'folded' : 
            style = '-' 
        else : 
            style = '--' 

        datafile = '%s/%s/%s'%(solvent,state,inFile) 
        try : 
            data = np.genfromtxt(datafile,skip_header=25) 
        except IOError : 
            print "No file found for %s/%s"%(solvent,state) 
            continue 

        data = data[equilTime:]

        data[:,0] /= 1000  ##ps -> ns 
        data[:,1] *= 10 #nm -> A
        
        #rc_file('%s/rc_files/%s'%(projectDir,rcFile) ) 
        
        ax.plot(data[:,0],data[:,1],color=color,linestyle=style) 
        print "%10s%10s%5.2f%5.2f"%(state, solvent,np.average(data[:,1]),np.std(data[:,1])) 

        binSize = 0.01
        hist, bin_edges = np.histogram(data[:,1],bins=np.arange(0,10,binSize),density=True) 

        ax2.plot(np.arange(0,10,binSize)[:-1], hist,color=color, linestyle=style) 
        

        #ax.set_xlim(0,xmax) 
        #ax.set_ylim(0.6,1.8) 
        #ax.set_title("%s %s"%(state,solvent)) 
    
fig.savefig("%s/rmsd_combined.png"%saveDir,format='png') 
plt.close(fig) 


fig2.savefig("%s/rmsd_histogram.png"%saveDir,format='png') 
plt.close(fig2) 
#print "%9s %9s plot complete"%(state,solvent)  





