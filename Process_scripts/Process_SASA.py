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
figRows=3

rc_file('rc_files/paper.rc') 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row',figsize=(3,6)) 
fig.subplots_adjust(left=0.15,right=0.95,top=0.95,bottom=0.10) 
fig.text(0.5,0.05, "Time", ha='center', va='center') 
fig.text(0.05,0.5, r"SASA", ha='center', va='center',rotation='vertical') 

binSize = 50 #ns
binSize *= 1000 #ns -> ps 
binSize /= 4 ##ps-> frames

index=0
xmax=50##Starting xlim
for solvent in ['water','tert','not_bound_sam'] : 
    ax = axarr[index%figRows]
    for state in ['folded','unfolded'] : 
        datafile = '%s/%s/dssp/totarea.xvg'%(solvent,state) 

        headlines = 0 
        with open(datafile) as f: 
            for line in f.readlines() : 
                if line.startswith('#') or line.startswith('@') : 
                    headlines += 1 
                else : 
                    break 

        data = np.genfromtxt(datafile,skip_header=headlines) 

        binNumber = len(data[:,0]) / binSize

        phobic =data[:,1]
        phillic=data[:,2]

        #Chop off data that does not fall into a bin
        while len(phobic) %binSize !=0 :
            phobic = phobic[:-1]
            phillic = phillic[:-1]
        assert len(phobic) == len(phillic)

        phobic = np.mean(phobic.reshape(-1,binSize),axis=1)
        phillic = np.mean(phillic.reshape(-1,binSize),axis=1)

        x = np.linspace(0,np.max(data[:,0]),binNumber)
        x /= 1000 

        if state == 'folded' : 
            c='b' 
        else :
            c='g'
        
        if data[:,0].max() > xmax : 
            xmax = data[:,0].max() 

        ax.plot(x,phillic,linestyle='--',color=c)
        ax.plot(x,phobic,linestyle='-',color=c)

    ax.set_title('%s'%(solvent.replace('_',' ')) ) 
    #ax.set_ylim(0,1)
    #ax.set_xlim(0,4.5)

    index +=1

#blu_dot = mlines.Line2D([],[], linestyle='None',color='b', marker='o',label="folded") 
#gre_dot = mlines.Line2D([],[], linestyle='None',color='g', marker='o',label="unfolded") 
#red_das = mlines.Line2D([],[],linestyle='--',color='r') 
#bla_das = mlines.Line2D([],[],linestyle='--',color='k') 

#fig.legend([blu_dot,gre_dot,red_das,bla_das],['folded','unfolded',r"R$_c$",r"2R$_c$"],numpoints=1,
#    loc='center', bbox_to_anchor=(0.9, 0.5),
#    fontsize='medium') 
fig.savefig('figures/sasa.png',format='png') 
