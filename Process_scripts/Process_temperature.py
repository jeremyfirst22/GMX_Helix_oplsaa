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

colorDict = {'water':'k','tert':'b','disordered_sam':'g'}

left, right = 0.13,0.95
bottom, top = 0.20,0.95

fig, ax = plt.subplots(1,1,figsize=(3,1.5))
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top)
fig.text((right-left)/2+left,0.00,           r"Time (ps)", ha='center', va='bottom')
fig.text(0.00,(top-bottom)/2+bottom,         r"Temperature (K)",ha='left',va='center',rotation='vertical')

index=0
xmax=50##Starting xlim
for solvent in ['water'] : #,'tert','disordered_sam'] : 
#    ax = axarr[index%figRows]
#    for state in ['folded','unfolded'] : 
        datafile = '%s/prep/Heating/temp/temperature.xvg'%solvent

        headlines = 0 
        with open(datafile) as f: 
            for line in f.readlines() : 
                if line.startswith('#') or line.startswith('@') : 
                    headlines += 1 
                else : 
                    break 
        try : 
            data = np.genfromtxt(datafile,skip_header=headlines) 
        except IOError : 
            print "No file found for %s"%(solvent) 
            continue

        ax.plot(np.arange(0,len(data[:,1])*0.2,0.2),data[:,1],color=colorDict[solvent]) #,label=state) 

ax.axvline(500,color='k',linestyle='--')
ax.axvline(600,color='k',linestyle='--')

ax.text(250,300,r"Heating",ha='center',va='center') 
ax.text(600+250,300,r"Cooling",ha='center',va='center') 

ax.set_xticks([0,500,600,1100]) 
fig.savefig('figures/temperature.png',format='png') 
