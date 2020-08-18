#!/usr/bin/env python 

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches
from matplotlib import rc_file
from scipy.signal import spline_filter
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
from scipy.signal import peak_widths
from sys import exit
from matplotlib.path import Path

rc_file('rc_files/paper.rc') 


##Green's theorem: https://stackoverflow.com/questions/22678990/how-can-i-calculate-the-area-within-a-contour-in-python-using-the-matplotlib
def area(vs) : 
    x = vs[:,0]
    y = vs[:,1]
    area=0.5*np.sum(y[:-1]*np.diff(x) - x[:-1]*np.diff(y))
    return np.abs(area)

colormap= plt.cm.get_cmap('jet')

 
for solvent in ['sam'] : #,'single_bound_sam','not_bound_sam'] : 
    for fold in ['folded'] :  #, 'unfolded'] : 
        left, right = 0.12, 0.92
        bottom, top =  0.08,0.98
        hspace, wspace = 0.25,0.15

        fig, ax = plt.subplots(1,1,figsize=(3.25,2))

        fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top)
        fig.subplots_adjust(wspace=wspace,hspace=hspace)
        fig.text(0.01,(top-bottom)/2+bottom,"Number of water molecules",rotation='vertical', ha='left', va='center') 

        if fold is 'folded' : 
            path=''
        elif fold is 'unfolded' : 
            path='../GMX_Helix_oplsaa_SAM_FROZEN/'
        print "%10s %10s"%(solvent, fold) 

        for index, selection in enumerate(['mainchain', 'sidechain', 'protein']) :
            datafile = "%s%s/%s/hbond/%s.xvg"%(path,solvent,fold,selection)
            headlines = 0 
            try : 
                with open(datafile) as f : 
                    for line in f.readlines() : 
                        if line.startswith('#') or line.startswith('@') : headlines += 1 
                        else : break 
                data = np.genfromtxt(datafile, skip_header=headlines) 
            except IOError : 
                print "ERROR: Data failed to import! %s"%datafile
                continue 

            datafile2 = "%s%s/%s/hbond/nearby_%s.xvg"%(path,solvent,fold,selection)
            headlines = 0 
            try : 
                with open(datafile2) as f : 
                    for line in f.readlines() : 
                        if line.startswith('#') or line.startswith('@') : headlines += 1 
                        else : break 
                data2 = np.genfromtxt(datafile2, skip_header=headlines) 
            except IOError : 
                print "ERROR: Data failed to import! %s"%datafile2
                continue 

            y1, y2 = data[:,1],data2[:,1]

            print "Num Nearby: %5.1f +/- %5.1f"%(np.average(y2), np.std(y2) ) 
            print "Num HBound: %5.1f +/- %5.1f"%(np.average(y1), np.std(y1) ) 

            width=0.35
            ax.bar(index, np.average(y1),color='gray',edgecolor='k',width=-width,align='edge') 
            ax.errorbar(index-0.5*width, np.average(y1),yerr=np.std(y1),color='k' ) 

            ax.bar(index, np.average(y2),color='b',edgecolor='k',width=width,align='edge') 
            ax.errorbar(index+0.5*width, np.average(y2),yerr=np.std(y2),color='k' ) 

        ax.set_xticks([0,1,2]) 
        ax.set_xticklabels(['Mainchain', 'Sidechain', 'Total protein']) 
        ax.set_ylim([0,125]) 


        blue = mpatches.Patch(color='b',label=r"Nearby water") 
        gray = mpatches.Patch(color='gray',label=r"Hydrogen bound water")
        ax.legend([blue,gray],["Nearby","Hydrogen bound"], loc=2)
        fig.savefig('figures/hbond_%s_%s.png'%(solvent,fold),format='png') 

