#!/usr/bin/env python

rcFile='paper.rc'
inFile ='rgyr/gyrate.xvg'
inFile2='rmsd/rmsd.xvg'
saveDir='figures'
outname='%s/combined_rgyr_v_rmsd.png'%saveDir

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm 
from matplotlib import rc_file
import glob
import os
import sys
from matplotlib import rc_file

figCols=2
figRows=3

rc_file("rc_files/%s"%rcFile) 

def Usage():
    print "Usage: %s <binSize>"%(sys.argv[0])

binSize=100
try : 
    binSize = int(sys.argv[1]) 
except : 
    Usage() 
    print "No bin size indicated, defaulting to %i"%binSize

binSize /= 4 #ps -> frames

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir) 

left, right = 0.08, 0.85
bottom, top =  0.1,0.95
hspace, wspace = 0.15,0.15

fig, axarr = plt.subplots(figRows, figCols, sharex='col',sharey='row', figsize=(5.0,3.8))
fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top)
fig.subplots_adjust(wspace=wspace,hspace=hspace)

fig.text((right-left)/2+left,0,           r"Radius of gyration R$_g$ ($\AA$)", ha='center', va='bottom')
fig.text(0,(top-bottom)/2+bottom,         r"RMSD ($\AA$)",ha='left',va='center',rotation='vertical')

fig.text((right-left)/4+left,top,            r"Starting with folded",ha='center',va='bottom')
fig.text(right-(right-left)/4,top,           r"Starting with unfolded",ha='center',va='bottom')
fig.text(right,top-(top-bottom)/6 ,          r"H$_2$O",ha='left',va='center',rotation=270)
fig.text(right,top-(top-bottom)/2,           r"2:1 H$_2$O:$t$-BuOH",ha='left',va='center',rotation=270)
fig.text(right,(top-bottom-hspace)/6+bottom, r"SAM surface",ha='left',va='center',rotation=270)

#legend

for row, solvent in enumerate(['water','tert','not_bound_sam']) : 
    for col,state in enumerate(['folded','unfolded']) : 
        if solvent == "water" : 
            equilTime = 600 
        else : 
            equilTime = 150 
        equilTime *= 1000 ## ns -> ps 
        equilTime /= 4    ## ps -> frames

        ax = axarr[row,col]

        datafile = "%s/%s/%s"%(solvent,state,inFile)  
        datafile2= "%s/%s/%s"%(solvent,state,inFile2) 
        
        try : 
            headlines = 0 
            with open(datafile) as f : 
                lines = f.readlines() 
                for line in lines : 
                    if line.startswith('#') or line.startswith('@') : 
                        headlines += 1 
                    else : 
                        break 
            data = np.genfromtxt(datafile,skip_header=headlines) 
        except IOError : 
            print "No file found for %s %s"%(state,solvent) 
            break  

        x = data[:,1] 
        x *= 10 ## nm -> AA
        avg1 = np.average(x[equilTime:]) 

        headlines=0 
        with open(datafile2) as f : 
            for line in f.readlines() : 
                if line.startswith('#') or line.startswith('@') : 
                    headlines += 1 
                else : 
                    break 

        data = np.genfromtxt(datafile2,skip_header=headlines,skip_footer=2)
        y = data[:,1]
        y *= 10 ## nm -> AA
        avg2 = np.average(y[equilTime:]) 

        while not len(x) == len(y) : 
            if len(x) > len(y) : 
                x = x[:-1]
            else : 
                y = y[:-1] 
        assert len(x) == len(y) 

        while len(x) % binSize !=0 : 
            x = x[:-1]
            y = y[:-1]
        assert len(x) % binSize == 0
        assert len(y) % binSize == 0

        xmin,xmax = 7.5, 17.5
        ymin,ymax = 1.50, 11 
#        xbins, ybins = np.arange(xmin,xmax,0.04),np.arange(ymin,ymax,0.04)
#        z, xbin, ybin = np.histogram2d(x ,y ,[xbins,ybins])
#        print np.max(z)
#        im = ax.pcolor(xbin,ybin,z.T,cmap='plasma',vmin=0,vmax=5000)
    
        ##Reshape into window averaged array
        xs = np.mean(x.reshape(-1,binSize),axis=1) 
        ys = np.mean(y.reshape(-1,binSize),axis=1) 
        
        #colors = cm.brg(np.linspace(0,1,len(ys)) ) 
        cm = plt.cm.get_cmap('brg') 
        z = np.arange(0,len(ys) ) 
        z = z*binSize*4/1000

        #counts,  xbins,  ybins  = np.histogram2d(x, y,range=((xmin,xmax),(ymin,ymax))) #,bins=(64,64))

        ax.plot(xs,ys,color='k',alpha=0.5,zorder=1,linewidth=1 ) 
        ax.set_xlim(xmin,xmax) 
        ax.set_ylim(ymin,ymax) 
        #ax.set_title('%s %s'%(state, solvent.replace('_',' '))) 


        
        sc = ax.scatter(xs ,ys,c=z,edgecolor='none',s=10,alpha=1,zorder=2,vmin=0,vmax=max(z)+1) #frames -> ns
        #sc = ax.scatter(np.linspace(0,1250,len(ys)),ys,c=z,edgecolor='none',s=2 , alpha=1, zorder=2, vmin=0, vmax=max(z)+1) 
        #print "Done with %s %s" %(solvent,state) 
        print "%15s%10s\t%5.2f\t%5.2f"%(solvent,state,avg1,avg2) 


cbar_ax = fig.add_axes([0.90, 0.15, 0.03, 0.75])
cbar = fig.colorbar(sc, cax=cbar_ax) 
#cbar.set_label('Time (ns)',rotation='vertical') 
fig.savefig(outname, format='png')
plt.close() 

