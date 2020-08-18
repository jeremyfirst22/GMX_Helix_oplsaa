#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob as glob 
from matplotlib import rc_file
from scipy.stats import linregress
from adjustText import adjust_text 
import os.path

probeStep = 0.01
probeMax  = 0.14 

rc_file('rc_files/paper.rc') 
colormap= plt.cm.get_cmap('viridis') 

figTA, axTA = plt.subplots(1,1)
Z = [[0,0],[0,0]]                                 # Apparently you can't set a colorbar without a mappable. This is hack
levels = np.arange(-probeStep/2,probeMax+probeStep/2, probeStep) #   around that. Taken from :
CS3 = axTA.contourf(Z, levels, cmap=colormap)     #   https://stackoverflow.com/questions/8342549/matplotlib-add-colorbar-to-a-sequence-of-line-plots
plt.close(figTA)


left, right = 0.15, 0.80
bottom, top = 0.15, 0.95

fig,ax = plt.subplots(1,1,figsize=(3.25,2)) 
fig.subplots_adjust(left=left, right=right, top=top,bottom=bottom) 
fig.text(0.01,(top-bottom)/2+bottom,"Number of waters", rotation='vertical', ha='left',va='center') 
fig.text((right-left)/2+left, 0.01, "$V_m$ (nm$^3$)", ha='center', va='bottom') 



for probe in np.arange(0.00, probeMax, probeStep) : 
    print "Probe = %f"%probe
    sizeAccum, volumeAccum = [], [] 
    for r in np.arange(0.05, 0.35, 0.02) : 
    
        file1 = "sam/folded/excluded_volume/size_%.2f.xvg"%r 
        try : 
            headlines = 0 
            with open(file1) as f : 
                for line in f.readlines() : 
                    if line.startswith('#') or line.startswith('@') : headlines += 1
                    else : break 
            data = np.genfromtxt(file1,skip_header=headlines) 
            size = np.average(data[:,1]) 
            sizeStd = np.std(data[:,1]) 
        except (IOError,IndexError) as e : 
            print "Data import error size %0.2f"%r
            print e
            continue 
    
        file2 = "sam/folded/excluded_volume/volume_p_%.2f_r_%.2f.xvg"%(probe,r) 
        try : 
            headlines = 0 
            with open(file2) as f : 
                for line in f.readlines() : 
                    if line.startswith('#') or line.startswith('@') : headlines += 1
                    else : break 
            data2 = np.genfromtxt(file2,skip_header=headlines) 
            volume = np.average(data2[:,1]) 
            volumeStd = np.std(data2[:,1]) 
        except (IOError,IndexError) as e : 
            print "Data import error volume %0.2f"%r
            print e
            continue 
    
        ax.scatter(volume,size,color=colormap(probe/probeMax)) 
        ax.errorbar(volume,size,xerr=volumeStd, yerr=sizeStd,color=colormap(probe/probeMax)) 

        volumeAccum.append(volume) 
        sizeAccum.append(size) 

    try : 
        slope, intercept, r_value, _, _ = linregress(volumeAccum, sizeAccum) 
    except : continue 
    
    xs = np.linspace(np.min(volumeAccum), np.max(volumeAccum),100) 
    ax.plot(xs, slope*xs+intercept, color=colormap(probe/probeMax)) 

cbar_ax = fig.add_axes([0.85, bottom, 0.03, top-bottom])
cbar = fig.colorbar(CS3, cax=cbar_ax,ticks=np.arange(0,probeMax,probeStep)) 
#cbar.ax.set_yticklabels(["0-%i ns"%(binSize/10),"%i-1250 ns"%(1250-(binSize/10))])

#ax.set_xlim([1.05,10.00]) 
#ax.set_ylim([0,150]) 
ax.vlines(5.7,80,100,color='k',linestyle='--',zorder=1000) 


fig.savefig('figures/volume_v_water.png',format='png') 




