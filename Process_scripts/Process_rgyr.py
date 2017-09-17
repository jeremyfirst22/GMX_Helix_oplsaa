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

for row,solvent in enumerate(['water','tert','sam']): 
    for col,state in enumerate(['folded','unfolded']) : 
        fig = plt.figure() 
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        fig.text(0.5,0.05, "Time (ns)", ha='center', va='center') 
        fig.text(0.05,0.5, r"R$_g$ (nm)",ha='center',va='center',rotation='vertical') 

        datafile = '%s/%s/%s'%(solvent,state,inFile) 
        try : 
            data = np.genfromtxt(datafile,skip_header=25) 
        except IOError : 
            print "No file found for %s/%s"%(solvent,state) 
            continue 
        data[:,0] /= 1000  ##ps -> ns 
        
        #rc_file('%s/rc_files/%s'%(projectDir,rcFile) ) 
        
        plt.plot(data[:,0],data[:,1]) 
        

        #ax.set_xlim(0,xmax) 
        ax.set_ylim(0.6,1.8) 
        ax.set_title("%s %s"%(state,solvent)) 
    
        fig.savefig("%s/rgyr_%s_%s.png"%(saveDir,state,solvent), format='png')
        plt.close(fig) 
        print "%9s %9s plot complete"%(state,solvent)  





