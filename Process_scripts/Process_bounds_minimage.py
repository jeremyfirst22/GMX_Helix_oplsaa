import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 

cutoff = 10 ##Electrostatic cutoff distance (Angstroms) 

figCols=1
figRows=3

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(left=0.1,right=0.8,top=0.9,bottom=0.1) 
fig.text(0.5,0.05, "Time (ns)", ha='center', va='center') 
fig.text(0.05,0.5, r"Distance to nearest image ($\AA$)", ha='center', va='center',rotation='vertical') 

index=0
xmax=50##Starting xlim
for solvent in ['sam','single_bound_sam','not_bound_sam'] : 
    ax = axarr[index%figRows]
    for state in ['folded','unfolded'] : 
        datafile = '%s/%s/minimage/mindist.xvg'%(solvent,state) 

        try : 
            headlines = 0 
            with open(datafile) as f : 
                lines = f.readlines() 
                for line in lines : 
                    if line.startswith('#') or line.startswith('@') : 
                        headlines += 1 
                    else : 
                        break 
            
            data = np.genfromtxt(datafile,skip_header=headlines) 
        except IOError : 
            print "No file found for %s %s"%(state,solvent)  
            continue
        if state == 'folded' : 
            c='b' 
        else :
            c='g'

        data[:,0] = data[:,0] / 1000 
        data[:,1] = data[:,1] * 10 # nm -> Angstroms
        
        #if data[:,0].max() > xmax : 
        xmax = data[:,0].max() 

        rvdw = np.full(len(data[:,1]),cutoff) 

        ax.scatter(data[:,0],data[:,1],s=0.1,color=c,label=state) 
        ax.plot(data[:,0],rvdw,'r--')
        ax.plot(data[:,0],2*rvdw,'k--') 
        ax.set_title('%s'%(solvent) ) 
        ax.set_xlim(0,xmax) 
        ax.set_ylim(0,60.0) 

    index +=1

blu_dot = mlines.Line2D([],[], linestyle='None',color='b', marker='o',label="folded") 
gre_dot = mlines.Line2D([],[], linestyle='None',color='g', marker='o',label="unfolded") 
red_das = mlines.Line2D([],[],linestyle='--',color='r') 
bla_das = mlines.Line2D([],[],linestyle='--',color='k') 

fig.legend([blu_dot,gre_dot,red_das,bla_das],['folded','unfolded',r"R$_c$",r"2R$_c$"],numpoints=1,
    loc='center', bbox_to_anchor=(0.9, 0.5),
    fontsize='medium') 
fig.savefig('figures/bounds_minimage.png',format='png') 
