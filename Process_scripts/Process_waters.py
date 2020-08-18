#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob as glob 
from matplotlib import rc_file
from scipy.stats import linregress
from adjustText import adjust_text 
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
import os.path

###
#  Define global variables
###
lowExpWater = 75   # Minimum number of waters (determined from Experiment) 
upExpWater  = 111  # Maximum number of waters (determined from Experiment)
rValues = np.arange(0.05, 0.35, 0.02) # Inclusion distances calculated with Gromacs
rValues *= 10                         # nm -> AA
splineKind = 'quadratic'              # Order of spline fit for interpolations
Vm    = 5.760      # Actual molecular volume (determined from Experiment) 
VmStd = 0.070      # Standard deviation of molecular volume (determined from Experiment) 

###
#  Given an input array and a target value, return the closest 
###
def find_nearest(array, value) : 
    array = np.asarray(array)
    idx = np.argmin(np.abs(array - value))
    return idx

###
#  Figure set up and labeling
###
rc_file('rc_files/paper.rc') 

left, right = 0.12, 0.98
bottom, top = 0.05, 0.98
hspace, wspace = 0.25, 0.00 

fig,axarr = plt.subplots(3,1,figsize=(3.25,6)) 
fig.subplots_adjust(left=left, right=right, top=top,bottom=bottom,hspace=hspace, wspace=wspace) 
fig.text(0.01,top-(top-bottom-hspace/2)/6,"Number of waters", rotation='vertical', ha='left',va='center') 
fig.text(0.01,top-(top-bottom)/2,"Radial distribution", rotation='vertical', ha='left',va='center') 
fig.text(0.01,bottom+(top-bottom-hspace/2)/6,"Volume ($\\times 10^3$ $\AA^3$)", rotation='vertical', ha='left',va='center') 
fig.text((right-left)/2+left, 0.0, "$d$ ($\AA$)", ha='center', va='bottom') 
fig.text((right-left)/2+left, 0.69, "$d$ ($\AA$)", ha='center', va='top') 
fig.text((right-left)/2+left, 0.355, "$d$ ($\AA$)", ha='center', va='top') 

ax1, ax2, ax3 = axarr 

ax1.text(0.05,0.95,r"\textsf{A}",va='top',ha='left',transform=ax1.transAxes,fontsize=12)
ax2.text(0.05,0.95,r"\textsf{B}",va='top',ha='left',transform=ax2.transAxes,fontsize=12)
ax3.text(0.05,0.95,r"\textsf{C}",va='top',ha='left',transform=ax3.transAxes,fontsize=12)

###
#  Subplot A: Number of waters as function of inclusion distance
###
ax1.axhline(lowExpWater, color='k', linestyle='--') 
ax1.axhline(upExpWater,color='k', linestyle='--') 
ax1.axhspan(lowExpWater, upExpWater,zorder=0, color='black',alpha=0.25,linewidth=0) 

sizeAccum = []
upLimAccum, lowLimAccum = [],[] 
for r in rValues : 
    file1 = "sam/folded/excluded_volume/size_%.2f.xvg"%(r/10.)  
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

    #ax1.scatter(r,size) 
    ax1.errorbar(r,size,yerr=sizeStd,color='blue',capsize=0) 

    sizeAccum.append(size) 
    upLimAccum.append(size+sizeStd) 
    lowLimAccum.append(size-sizeStd) 

sizeAccum = np.array(sizeAccum,dtype=float) 
spl = interp1d(rValues, sizeAccum, kind=splineKind) 

xs = np.linspace(np.min(rValues), np.max(rValues), 1000) 
ax1.plot(xs, spl(xs),'blue' ) 

upLim = interp1d(rValues, upLimAccum, kind=splineKind)
lowLim = interp1d(rValues, lowLimAccum, kind=splineKind) 
ax1.fill_between(xs, upLim(xs), lowLim(xs), color='blue', alpha=0.5,linewidth=0) 

minInclusDist, maxInclusDist = xs[find_nearest(spl(xs), lowExpWater)], xs[find_nearest(spl(xs),upExpWater)]
    ##Minimum and maximum inclusion distances that correspond with experimentally determined waters
ax1.axvline(minInclusDist, color='k', linestyle=':')  
ax1.axvline(maxInclusDist, color='k', linestyle=':') 
ax1.set_xlim([0.45,3.35]) 

print minInclusDist, maxInclusDist

###
#  Subplot B: Radial distribution functions
###
#fileList = ["sam/folded/rdf/protein.xvg", "sam/folded/rdf/lys.xvg", "sam/folded/rdf/leu.xvg", "sam/folded/rdf/back.xvg", "sam/folded/rdf/lys_wat.xvg"] #, "sam/folded/rdf/backNH_wat.xvg"] 
#fileList = [ "sam/folded/rdf/lys_wat.xvg", "sam/folded/rdf/lys_test.xvg"] 
#fileList = ["sam/folded/rdf/lys_wat.xvg", "sam/folded/rdf/leu_wat.xvg", "sam/folded/rdf/backNH_wat.xvg", "sam/folded/rdf/backCO_wat.xvg", "sam/folded/rdf/protein_test.xvg"] 
fileList = ["sam/folded/rdf/protein.xvg"] 
for file2 in fileList :     
    if "protein" in file2 : 
        name = "Protein" 
        color = 'k' 
    elif "lys" in file2 : 
        name = "Lys" 
        color = 'r' 
    elif "leu" in file2 : 
        name = "Leu" 
        color = 'c' 
    elif "back" in file2 : 
        name = "Backbone" 
        color = 'g' 

    headlines = 0 
    with open(file2) as f: 
        for line in f.readlines() : 
            if line.startswith('#') or line.startswith('@') : 
                headlines += 1 
            else : break 
    data2 = np.genfromtxt(file2, skip_header=headlines) 
    data2[:,0] *= 10 ## nm -> AA

    if "lys_wat" in file2 : data2[:,0] -= 1.01
    if "leu_wat" in file2 : data2[:,0] -= 1.08
    if "backNH_wat" in file2 : data2[:,0] -= 1.01
    if "protein" in file2 : data2[:,0] -= 1.05
    
    ax2.plot(data2[:,0], data2[:,1], color=color, label=name) 

#ax2.legend(loc=1) 
ax2.axvline(minInclusDist, color='k', linestyle=':')  
ax2.axvline(maxInclusDist, color='k', linestyle=':') 

ax2.set_xlim([0.5, 8.00]) 
#ax2.set_xlim([0.45,3.35]) 

###
#  Subplot C: Volume as a function of probe radius 
###
probeMax, probeStep = 0.14, 0.01 
colorList = ['b', 'g', 'r'] 
for c, probe in enumerate([0.00, 0.03, 0.14]) : 
    volumeAccum = [] 
    upLimAccum, lowLimAccum = [], [] 
    if probe == 0 : label = "VDW volume"
    elif probe == 0.14 : label = "SA volume"
    else  : label = "Probe of %.1f $\AA$"%(probe*10)

    for i,r in enumerate(rValues) : 
        file2 = "sam/folded/excluded_volume/volume_p_%.2f_r_%.2f.xvg"%(probe,r/10.)
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

        ax3.errorbar(r,volume,yerr=volumeStd, color=colorList[c],capsize=0) 

        volumeAccum.append(volume)
        upLimAccum.append(volume+volumeStd)
        lowLimAccum.append(volume-volumeStd)

    volumeAccum = np.array(volumeAccum,dtype=float) 
    spl = interp1d(rValues, volumeAccum, kind=splineKind)

    xs = np.arange(np.min(rValues), np.max(rValues), 0.0001) 
    ax3.plot(xs, spl(xs),color=colorList[c],label=label ) 

    upLim = interp1d(rValues, upLimAccum, kind=splineKind) 
    lowLim = interp1d(rValues, lowLimAccum, kind=splineKind) 
    ax3.fill_between(xs, upLim(xs), lowLim(xs),color=colorList[c],alpha=0.5, linewidth=0) 

ax3.axvline(minInclusDist, color='k', linestyle=':')  
ax3.axvline(maxInclusDist, color='k', linestyle=':') 

ax3.axhline(Vm-VmStd,color='k', linestyle='--') 
ax3.axhline(Vm+VmStd,color='k', linestyle='--') 
ax3.axhspan(Vm+VmStd, Vm-VmStd, zorder=0, color='black', alpha=0.25, linewidth=0) 
ax3.legend(loc=(0.12,0.65))#,fontsize=6) 
ax3.set_xlim([0.45,3.35]) 
#ax3.set_ylim([1.5, 14.5]) 

fig.savefig('figures/water.png',format='png') 

