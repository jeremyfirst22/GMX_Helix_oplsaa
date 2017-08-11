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

def Usage():
    print "Usage: %s <Data directory>"%(sys.argv[0])

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir) 

curdir = os.path.abspath(os.getcwd())

sol2inFileList={
'tert':['tba_tba.xvg','tba_wat.xvg','wat_wat.xvg','leu_tba.xvg','leu_wat.xvg','lys_tba.xvg','lys_wat.xvg'],
'water':['wat_wat.xvg','leu_wat.xvg','lys_wat.xvg'],
'sam':['wat_wat.xvg','leu_wat.xvg','lys_wat.xvg','sam_leu.xvg','sam_lys.xvg']
} 

for sol in 'tert', 'water', 'sam' : 
    for inFile in sol2inFileList[sol] : 
        title=os.path.basename(inFile).split('.')[0]
        
        outname = "%s/rdf_%s_%s.pdf"%(saveDir,sol,title)
        outname = os.path.join(curdir,outname)
        
        for state in 'folded', 'unfolded' : 
            datafile = "%s_%s/rdf/%s"%(state,sol,inFile) 
            try : 
                data = np.genfromtxt(datafile,skip_header=16) 
            except IOError :
                print "Error: %s not found."%datafile
                break
            except ValueError : 
                print "Error: skip_header lines incorrect for %s"%datafile
                sys.exit()  
            x = data[:,0] 
            y = data[:,1] 
        
            plt.plot(x ,y) 
        else : 
            plt.title(title) 
            plt.xlabel(r"r (nm)") 
            plt.ylabel(r"p(r)") 
            plt.xlim([0,2]) 
            plt.gca().set_ylim(bottom=0) 

            blu_lin = mlines.Line2D([],[],linestyle='-',color='b',label="folded") 
            gre_lin = mlines.Line2D([],[],linestyle='-',color='g',label="unfolded") 
            plt.legend(handles=[gre_lin,blu_lin],loc=1,fontsize='medium') 
        
            plt.savefig(outname, format='pdf')
            print "Success: %s_%s plot saved to %s"%(state,sol,saveDir) 
        plt.close()

