import numpy as np 
import matplotlib.pyplot as plt 
import sys 
from matplotlib import rc_file

ymax = {
'tba_tba':3,
'tba_wat':1.4,
'wat_wat':6 
}

rc_file('rc_files/paper.rc') 


figRows, figCols = 3,1
left, right = 0.15,0.95
bottom, top= 0.05,0.98
fig, axarr = plt.subplots(figRows, figCols, sharex='all',sharey='row', figsize=(3.0,6.0))
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top)
fig.subplots_adjust(wspace=0,hspace=0.1)

fig.text((right-left)/2+left,0,           r"g(r)", ha='center', va='bottom')
fig.text(0,(top-bottom)/2+bottom,         r"r ($\AA$)",ha='left',va='center',rotation='vertical')

for index,item in enumerate(['tba_tba', 'tba_wat', 'wat_wat']) : 
    ax = axarr[index]
    ax.text(0.95,0.95,item.replace('_','-').replace('tba','TBA').replace('wat','H$_2$O'),transform=ax.transAxes,ha='right',va='top')

    fileName = "rdfs_for_tert_validation/%s.xvg"%item
    headlines = 0 
    with open(fileName) as f : 
        for line in f.readlines() : 
            if line.startswith('#') or line.startswith('@') : 
                headlines += 1 
            else : break 
    data = np.genfromtxt(fileName,skip_header=headlines)

    data[:,0] *= 10
    y = np.full_like(data[:,0],1)

    ax.plot(data[:,0],data[:,1],color='m',linewidth=2)#,s=0.1 
    ax.plot(data[:,0],y,color='k') 
    
    ax.set_ylim(0,ymax[item]) 
    ax.set_xlim(0,20)
    
plt.savefig("figures/rdfs.png",format='png') 

sys.exit() 
