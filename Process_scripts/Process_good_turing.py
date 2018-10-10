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
fig.text(0.5,0.05, "RMSD ($\AA$)", ha='center', va='center') 
fig.text(0.05,0.5, r"Prob of unobserved", ha='center', va='center',rotation='vertical') 

index=0
xmax=50##Starting xlim
for solvent in ['water','tert','sam'] : 
    ax = axarr[index%figRows]
    for state in ['folded','unfolded'] : 
        datafile = '%s/%s/good-turing/good_turing.rmsd.P_unobserved_vs_RMSD.dat'%(solvent,state) 

        try : 
            
            data = np.genfromtxt(datafile,skip_header=1) 
        except IOError : 
            print "No file found for %s %s"%(state,solvent)  
            continue
        if state == 'folded' : 
            c='b' 
        else :
            c='g'
        
        if data[:,0].max() > xmax : 
            xmax = data[:,0].max() 

        ax.errorbar(data[:,0],data[:,1],yerr=data[:,2],color=c,label=state) 
        ax.set_title('%s'%(solvent) ) 

    index +=1

#blu_dot = mlines.Line2D([],[], linestyle='None',color='b', marker='o',label="folded") 
#gre_dot = mlines.Line2D([],[], linestyle='None',color='g', marker='o',label="unfolded") 
#red_das = mlines.Line2D([],[],linestyle='--',color='r') 
#bla_das = mlines.Line2D([],[],linestyle='--',color='k') 

#fig.legend([blu_dot,gre_dot,red_das,bla_das],['folded','unfolded',r"R$_c$",r"2R$_c$"],numpoints=1,
#    loc='center', bbox_to_anchor=(0.9, 0.5),
#    fontsize='medium') 
fig.savefig('figures/good_turing.png',format='png') 
