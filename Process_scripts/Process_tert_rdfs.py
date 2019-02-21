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


figRows, figCols = 1,1
left, right = 0.10,0.95
bottom, top= 0.15,0.95
fig, ax    = plt.subplots(figRows, figCols, sharex='all',sharey='row', figsize=(3.25,2.00))
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top)
fig.subplots_adjust(wspace=0,hspace=0.1)

fig.text((right-left)/2+left,0.01,           r"g(r)", ha='center', va='bottom')
fig.text(0.01,(top-bottom)/2+bottom,         r"r ($\AA$)",ha='left',va='center',rotation='vertical')

colors = {
        'tba_tba':'b',
        'tba_wat':'g',
        'wat_wat':'r' 
        }

for index,item in enumerate(['wat_wat', 'tba_wat', 'tba_tba']) : 
#    ax = axarr[index]
    #ax.text(0.95,0.95,item.replace('_','-').replace('tba','TBA').replace('wat','H$_2$O'),transform=ax.transAxes,ha='right',va='top')

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

    ax.plot(data[:,0],data[:,1],linewidth=2,
            color=colors[item],
            label=item.replace('_','-').replace('tba','TBA').replace('wat','H$_2$O'))#,s=0.1 
    ax.plot(data[:,0],y,color='k') 
    
    ax.set_ylim(0,6) 
    ax.set_xlim(0,20)

ax.legend(loc=1) 
    
plt.savefig("figures/rdfs.png",format='png') 

sys.exit() 
