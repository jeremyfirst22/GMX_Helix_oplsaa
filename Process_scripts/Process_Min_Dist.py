import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 

figCols=1
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0) 
fig.text(0.5,0.04, "Time (ns)", ha='center', va='center') 
fig.text(0.08,0.5, r"Distance to nearest image ($\AA$)", ha='center', va='center',rotation='vertical') 

index=0
for solvent in ['water','tert','sam'] : 
    ax = axarr[index%figRows]
    for state in ['folded','unfolded'] : 
        datafile = '%s_%s/minimage/mindist.xvg'%(state,solvent) 

        try : 
            data = np.genfromtxt(datafile,skip_header=27) 
        except IOError : 
            print "No file found for %s %s"%(state,solvent)  
            continue
        if state == 'folded' : 
            c='b' 
        else :
            c='g'

        data[:,0] = data[:,0] / 1000 
        data[:,1] = data[:,1] * 10 # nm -> Angstroms

        ax.scatter(data[:,0],data[:,1],s=0.1,color=c,label=state) 
        ax.set_title('%s'%(solvent) ) 
        ax.set_xlim(0,50) 
        ax.set_ylim(10,50.0) 

    index +=1

#blu_dot = mlines.Line2D([],[], linestyle='None',color='b', marker='o',label="folded") 
#gre_dot = mlines.Line2D([],[], linestyle='None',color='g', marker='o',label="unfolded") 

#plt.legend(handles=[gre_dot,blu_dot],numpoints=1,loc='center left', bbox_to_anchor=(1.0, 1.75))
fig.savefig('figures/minimage.png',format='png') 
