import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
from matplotlib import rc_file
import matplotlib.lines as mlines 

cutoff = 10 ##Electrostatic cutoff distance (Angstroms) 

figCols=1
figRows=1

rc_file('rc_files/paper.rc') 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

colorDict = {'water':'k','tert':'b','not_bound_sam':'g'}

left, right = 0.13,0.95
bottom, top = 0.13,0.95

fig, ax = plt.subplots(1,1)
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top)
fig.text((right-left)/2+left,0.00,           r"RMSD ($\AA$)", ha='center', va='bottom')
fig.text(0.00,(top-bottom)/2+bottom,         r"Probability of unobserved structure",ha='left',va='center',rotation='vertical')

index=0
xmax=50##Starting xlim
for solvent in ['water','tert','not_bound_sam'] : 
#    ax = axarr[index%figRows]
#    for state in ['folded','unfolded'] : 
        datafile = '%s/good_turing/good_turing.rmsd.P_unobserved_vs_RMSD.dat'%(solvent) 

        try : 
            
            data = np.genfromtxt(datafile,skip_header=1) 
        except IOError : 
            print "No file found for %s"%(solvent) 
            continue
        #if state == 'folded' : 
        #c='b' 
        #else :
        #    c='g'
        
        if data[:,0].max() > xmax : 
            xmax = data[:,0].max() 

        try : 
            ax.errorbar(data[:,0],data[:,1],yerr=data[:,2],color=colorDict[solvent]) #,label=state) 
        except IndexError : 
            ax.scatter(data[:,0],data[:,1],s=0.5,color=colorDict[solvent]) #,label=state) 
            ax.plot(data[:,0],data[:,1],color=colorDict[solvent]) #,label=state)    
    #ax.set_title('%s'%(solvent.replace('_',' ')) ) 
ax.set_ylim(0,1)
ax.set_xlim(0,4.5)

#index +=1

#blu_dot = mlines.Line2D([],[], linestyle='None',color='b', marker='o',label="folded") 
#gre_dot = mlines.Line2D([],[], linestyle='None',color='g', marker='o',label="unfolded") 
#red_das = mlines.Line2D([],[],linestyle='--',color='r') 
#bla_das = mlines.Line2D([],[],linestyle='--',color='k') 

#fig.legend([blu_dot,gre_dot,red_das,bla_das],['folded','unfolded',r"R$_c$",r"2R$_c$"],numpoints=1,
#    loc='center', bbox_to_anchor=(0.9, 0.5),
#    fontsize='medium') 
fig.savefig('figures/good_turing.png',format='png') 
