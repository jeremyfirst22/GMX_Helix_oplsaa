#!/usr/bin/env python

rcFile='paper.rc'
saveDir='figures'

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm 
import glob
import os
import sys
from matplotlib import rc_file
import matplotlib.lines as mlines 

rcFile = 'rc_files/presentation.rc' 

def Usage():
    print "Usage: %s <Data directory>"%(sys.argv[0])

rc_file(rcFile) 

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir) 

curdir = os.path.abspath(os.getcwd())

sol2inFileList={
'tert':['leu_tba.xvg','leu_wat.xvg','lys_tba.xvg','lys_wat.xvg'],
'water':['leu_wat.xvg','lys_wat.xvg'],
'sam':['leu_wat.xvg','lys_wat.xvg','sam_leu.xvg','sam_lys.xvg']
} 

for sol in 'tert', 'water', 'sam' : 
    for inFile in sol2inFileList[sol] : 
        title=os.path.basename(inFile).split('.')[0]
        
        outname = "%s/rdf_%s_%s.png"%(saveDir,sol,title)
        outname = os.path.join(curdir,outname)
        
        for state in 'folded', 'unfolded' : 
            if state == 'folded' : 
                color = 'b' 
            else : 
                color = 'g' 
            datafile = "%s/%s/rdf/%s"%(sol,state,inFile) 
            try : 
                headlines = 0 
                with open(datafile) as f : 
                    for line in f.readlines() : 
                        if line.startswith('#') or line.startswith('@') : 
                            headlines += 1
                        else : 
                            break 
                data = np.genfromtxt(datafile,skip_header=headlines) 
            except IOError :
                print "Error: %s not found."%datafile
                break
            except ValueError : 
                print "Error: skip_header lines incorrect for %s"%datafile
                sys.exit()  
            x = data[:,0] 
            y = data[:,1] 
        
            plt.plot(x ,y,color=color) 
        else : 
            plt.title(title) 
            plt.xlabel(r"r (nm)") 
            plt.ylabel(r"p(r)") 
            plt.xlim([0,2]) 
            plt.gca().set_ylim(bottom=0) 

            blu_lin = mlines.Line2D([],[],linestyle='-',color='b',label="folded") 
            gre_lin = mlines.Line2D([],[],linestyle='-',color='g',label="unfolded") 
            plt.legend(handles=[gre_lin,blu_lin],loc=1,fontsize='medium') 
        
            plt.savefig(outname, format='png')
            print "Success: %s_%s plot saved to %s"%(state,sol,saveDir) 
        plt.close()

