#!/usr/bin/env python 

import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc_file
from scipy.signal import spline_filter
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
from scipy.signal import peak_widths
import matplotlib.gridspec as gridspec
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
        left, right = 0.15, 0.92
        bottom, top =  0.08,0.98
        hspace, wspace = 0.25,0.15

        cbSpace = 0.10

        fig = plt.figure(figsize=(3.25,4))
        axarr = gridspec.GridSpec(2,1,height_ratios=[1, 2.00],wspace=wspace,hspace=hspace,left=left,right=right,bottom=bottom,top=top) 

        #fig, axarr = plt.subplots(2,1,figsize=(3.25,  4.75), sharex='none', sharey='none') 

        #fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top)
        #fig.subplots_adjust(wspace=wspace,hspace=hspace)
        
        fig.text(0.04,0.35, r"$x$ ($\AA$)",ha='left',va='center',rotation='vertical')
        fig.text((right-left-cbSpace)/2+left,0.01  , r"$y$          ($\AA$)",ha='center',va='bottom')                   
        fig.text((right-left)/2+left,0.64  , r"$z$          ($\AA$)",ha='center',va='bottom')                   
        fig.text(0.04,0.84,    r"$\rho_{\rm{number}}$",ha='left',va='center',rotation='vertical')
        
        ax1 = plt.subplot(axarr[0])  
        ax2 = plt.subplot(axarr[1])  

        ax1.text(0.035,0.90,r"\textsf{A}",va='top',ha='left',transform=ax1.transAxes,fontsize=12)
        ax2.text(0.05,0.95,r"\textsf{B}",va='top',ha='left',transform=ax2.transAxes,fontsize=12,color='w')

        box = ax2.get_position()
        box.x1 -= cbSpace
        ax2.set_position(box) 
        ax2.set_aspect('equal',anchor='SW') 
        #ax1.set_aspect('equal') 

        if fold is 'folded' : 
            path=''
        elif fold is 'unfolded' : 
            path='../GMX_Helix_oplsaa_SAM_FROZEN/'
        print "%10s %10s"%(solvent, fold), 

        peakHeight = 0.01
        datafile = "%s%s/%s/density/density.xvg"%(path,solvent,fold)
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

        x, y = data[:,0],data[:,1]
        x *= 10  ##nm -> A

        ax1.plot(x,y,color='k') 

        peaks, properties = find_peaks(y,height=8,prominence=1) 

        results = peak_widths(y, peaks, rel_height=1-peakHeight)  

        maxWidth = np.argmax([results[0]]) 
        xMin = x[int(np.round(results[2][maxWidth]))]  
        xMax = x[int(np.round(results[3][maxWidth]))]  

        #plt.fill_between(x,results[1][maxWidth],y, y > results[1][maxWidth] ) 

        ax1.hlines(results[1][maxWidth],xMin, xMax, color='k', linestyle=':') 
        #plt.hlines(results_half[maxWidth][0], x[results[maxWidth][1], x[results[maxWidth][0]) 

#        plt.hlines(*results_half[1:], color="C2",transform=plt.gca().get_yaxis_transform()) 
#        ys = results_half[1:][0]
#        xmins = results_half[1:][1]
#        xmaxs = results_half[1:][2]
#        for i in range(len(ys)) : 
#            plt.hlines(ys[i],x[xmins[i]],x[xmaxs[i]]) 

        print "Height = %.3f A"%(xMax-xMin) 
        ax1.text(0.95,0.95,r"Height = %.1f $\AA$"%(xMax-xMin),transform=ax1.transAxes,ha='right',va='top') 

        peakHeight = 0.25 

        datafile = "%s%s/%s/density/densmap.dat"%(path,solvent,fold)
        with open(datafile) as f :
            lines = f.readlines()
            for line in lines[:2] :
                if 'x-range:' in line :
                    xMin = float(line.split()[2])
                    xMax = float(line.split()[3])
                if 'y-range:' in line :
                    yMin = float(line.split()[2])
                    yMax = float(line.split()[3])

        data = np.genfromtxt(datafile,skip_header=2,dtype=float)

        #data -= np.max(data)    ##Set max dG (not sampled region) to zero. dG are relative to that.

        #data[data == 0] = 'nan' ##Replace dG = 0 with 'nan' so it doesn't show up on heat map.

        x = np.linspace(xMin,xMax,np.shape(data)[1])
        x *= 10 ## nm -> A
        y = np.linspace(yMin,yMax,np.shape(data)[0])
        y *= 10 ## nm -> A

        data /= 100 ## counts -> 10^3 counts 
        ##
        #   z-axis unit conversion not needed. Units are actually in absolute counts
        #data = data / 100 ## count/nm^2 -> count/A^2
        ##

        vmin, vmax = 0,np.max(data) 
        im = ax2.pcolor(x,y, data, cmap=colormap,vmin=vmin, vmax=vmax) 

        X, Y = np.meshgrid(x,y) 
        #levels = [100, 500, 1000, 2000, 4000] 
        levels = [np.max(data)*peakHeight]
        cs = ax2.contour(X,Y,data,levels=levels,colors='white') 
        #ax2.clabel(cs, inline=1, fontsize=10,colors='white')

        X, Y = X.flatten(), Y.flatten() 
        points = np.vstack((X,Y)).T

        for i in range(len(levels)):
            contour = cs.collections[i]
            vs = contour.get_paths()[0].vertices
            p = Path(vs) 
            grid = p.contains_points(points) 
            mask = grid.reshape(np.size(y), np.size(x) ) 
            #print np.size(mask) 
            #grid = p.contains_points(
            # Compute area enclosed by vertices.
            a = area(vs)
            #print "r = " + str(levels[i]) + ": a =" + str(a)
            print "Area = %.3f .  Num enclosed pixels: %i"%(a,np.sum(mask)) 
            pixArea = np.average(np.diff(x)) * np.average(np.diff(y)) 
            print "\tArea from pixels = %.3f A^2"%(np.sum(mask) * pixArea) 
        ax2.text(0.95,0.95,r"Area = %.1f $\AA^2$"%(np.sum(mask) * pixArea),transform=ax2.transAxes,ha='right',va='top',color='white') 
    
        #for i in np.arange(0,np.shape(mask)[0],5)  : 
        #    for j in np.arange(0,np.shape(mask)[1],5) : 
        #        print i, j 
        #        if mask[i,j] : ax2.scatter(x[j], y[i], s = 1.0,color='white') 

        #figTA, axTA = plt.subplots(1,1)
        #Z = [[0,0],[0,0]]                                 # Apparently you can't set a colorbar without a mappable. This is hack
        #levels = np.linspace(vmin,vmax,numLevels) #   around that. Taken from :
        #CS3 = axTA.contourf(Z, levels, cmap=colormap)     #   https://stackoverflow.com/questions/8342549/matplotlib-add-colorbar-to-a-sequence-of-line-plots
        #plt.close(figTA)

        leftbar = right - 2*cbSpace/3
        widthbar = cbSpace*1/3

        cbar_ax = fig.add_axes([leftbar, bottom, widthbar, ax2.get_position().y1 - bottom]) 
        cbar = fig.colorbar(im, cax=cbar_ax,orientation='vertical')
        fig.text(leftbar, (ax2.get_position().y1 - bottom) / 2 + bottom, r"Counts $\times 10^3$", rotation='vertical',ha='right',va='center') 

        ax1.set_xlim([10,45]) 
        ax2.set_xlim([10,50]) 
        ax2.set_ylim([15,55]) 

        #ax2.hlines(35, 22.5, 38.5, color='w', linestyle=':') 
        #ax2.axvline(22.5, color='w', linestyle=':') 
        #ax2.axvline(38.5, color='w', linestyle=':') 

        fig.savefig("figures/dimensions_%s_%s.png"%(solvent, fold), format='png') 
