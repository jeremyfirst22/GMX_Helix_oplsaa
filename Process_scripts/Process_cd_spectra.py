#!/usr/bin/env python 

projectDir='Users/jeremyfirst/GMX_Helix_oplsaa'
#rcFile='paper.rc'
rcFile='paper.rc'
saveDir='figures'

import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc_file
import glob as glob 
import os 
import sys 

rc_file('rc_files/%s'%rcFile) 
colormap= plt.cm.get_cmap('viridis')
figRows, figCols = 3,2

binSize= 125.0 ##ns 
skip = 1     ##Skip every n frames 

binSize*=10  ##ns->framess
binSize = int(binSize) 

if not binSize%skip == 0 : 
    print "Error: Skip/binSize wrong size. Skip does not fit into binSize" 
    sys.exit() 

top=os.getcwd()

colorDict = {'water':'k','tert':'b','not_bound_sam':'g'}

xmin,xmax, xstep = 190,250,1.0
#bandwidth = 10.5 

indexToLetter = {
        1:'A',
        2:'B',
        3:'C',
        4:'D',
        5:'E',
        6:'F'
        }


def gauss (a,b,c,x) :
    return a*np.exp(-(x-b)**2/(c**2))

for bandwidth in [10.0] : #np.arange(9,13,0.5,dtype=float) : 
    sys.stdout.write("\nBandwidth = %3.1f:\n"%bandwidth)  

    left, right = 0.15,0.98
    bottom, top = 0.13,0.95
    wspace = 0 

    fig, axarr = plt.subplots(1,2,figsize=(3.2,2.5),sharex='all',sharey='all') 
    fig.subplots_adjust(left=left, bottom=bottom,right=right,top=top,wspace=wspace)
    fig.text((right-left)/2+left,0.00,           r"Wavelength $\lambda $(nm)", ha='center', va='bottom')
    fig.text(0.01,(top-bottom)/2+bottom,         r"Ellipticity $\thetaup$ (10$^{3}$ deg cm$^2$ dmol$^{-1}$)",ha='left',va='center',rotation='vertical')

    for combCol, param in enumerate([ 'hirst' ,'woody' ]) : 
        sys.stdout.write("\nUsing %6s parameters:\n"%param)  

        ax = axarr[combCol]

        figTA, axTA = plt.subplots(1,1) 
        Z = [[0,0],[0,0]]                                 # Apparently you can't set a colorbar without a mappable. This is hack 
        step = 25 ##size of each step in color bar
        levels = np.arange(0,1250 + binSize/ 10 ,binSize/ 10) #   around that. Taken from :
        CS3 = axTA.contourf(Z, levels, cmap=colormap)     #   https://stackoverflow.com/questions/8342549/matplotlib-add-colorbar-to-a-sequence-of-line-plots
        plt.close(figTA) 


        left, right = 0.10, 0.75
        bottom, top =  0.1,0.95
        hspace, wspace = 0.15,0.1

        fig2, axarr2 = plt.subplots(figRows, figCols, sharex='all',sharey='all', figsize=(5.0,3.8))
        fig2.subplots_adjust(left=left, bottom=bottom,right=right,top=top)
        fig2.subplots_adjust(wspace=wspace,hspace=hspace)

        fig2.text((right-left)/2+left,0.03,           r"Wavelength $\lambda $(nm)", ha='center', va='center')
        fig2.text(0.02,(top-bottom)/2+bottom,         r"Ellipticity $\thetaup$ (10$^{3}$ deg cm$^2$ dmol$^{-1}$)",ha='center',va='center',rotation='vertical')
        
        fig2.text((right-left)/4+left,top,            r"Starting with folded",ha='center',va='bottom')
        fig2.text(right-(right-left)/4,top,           r"Starting with unfolded",ha='center',va='bottom')
        fig2.text(right,top-(top-bottom)/6 ,          r"H$_2$O",ha='left',va='center',rotation=270)
        fig2.text(right,top-(top-bottom)/2,           r"2:1 H$_2$O:$t$-BuOH",ha='left',va='center',rotation=270)
        fig2.text(right,(top-bottom-hspace)/6+bottom, r"SAM surface",ha='left',va='center',rotation=270)

        index = 0 
        for row,solvent in enumerate(['water', 'tert', 'not_bound_sam']) : 
            if solvent == 'water' :
                equilTime = 600 ##ns
            else :
                equilTime = 150 ##ns
            equilTime *= 10 # ns -> frames, frames printed every 100 ps

            color = colorDict[solvent] 
    
            xs = np.arange(xmin,xmax,xstep,dtype='float') 
            avgData = np.zeros_like(xs) 
    
            for col,state in enumerate(['folded', 'unfolded']) : 
                index += 1 
                ax2 = axarr2[row,col]

                datafiles = glob.glob('%s/%s/cd_spectra/output-%s-%s-%s/frame*.cdl'%(solvent,state,solvent.replace('_','-'),state,param) )
                numBins = int(len(datafiles) / binSize ) 
                
                numFiles = len(datafiles[equilTime:]) /skip    ##number of line spectra in average spectra
                for j in range(0,numBins) :    ## For each bin
                    sys.stdout.write('\r')
                    sys.stdout.write("Processing %5i of %5i %6s %8s"%(j+1,numBins,solvent,state) ) 
                    sys.stdout.flush()     

                    bindata = np.zeros_like(xs) 
                    for i in range(j*binSize,(j+1)*binSize,skip) : ##For each spectra in each bin
                        datafile = "%s/%s/cd_spectra/output-%s-%s-%s/frame%i.cdl"%(solvent,state,solvent.replace('_','-'),state,param,i)
                        data = np.genfromtxt(datafile) 
    
                        for trans in data :
                            a, b, c = trans[1], trans[0], bandwidth
                            spec = gauss(a,b,c,xs)
        
                            #plt.scatter(trans[0],trans[1])
                            #plt.plot(xs,spec)
        
                            for k,x in enumerate(xs) :
                                bindata[k] += spec[k]
                                if i > equilTime : 
                                    avgData[k] += spec[k]
                    bindata /= binSize / skip

                    ##This is an unit conversion from DBR (Rotational strength) to ellipticity. 
                    ##The 0.248 factor is from Schellman. Chem. Reviews. 1974 pg. 326. 
                    ##The 0.0082 factor is fit to match my line spectra with output from DichroCalc.
                    ## I'm not sure why I'm off by this factor, but is likely conversion between 
                    ##   delta epsilon (extinction) and molar ellipticity
                    bindata = (bindata * 0.248 * xs)/0.0082 
                    bindata /= 1000 

                    sc = ax2.plot(xs,bindata)#,label=label)


                #ax2.set_title("%s %s"%(solvent.replace('_',' '),state)) 
                ax2.set_ylim([-35,40]) 
                ax2.set_xlim([190,250]) 

                ax2.text(0.97,0.95,r"\textsf{%c}"%indexToLetter[index],va='top',ha='right',transform=ax2.transAxes,fontsize=12)

                colors = [colormap(i) for i in np.linspace(0,1,len(ax2.lines))] 
                for i,j in enumerate(ax2.lines) : 
                    j.set_color(colors[i]) 
       #         fig2.legend(loc=1) 

                sys.stdout.write('\n') 

            avgData /= numFiles * 2 ##two trajectories
            ##DBR -> mdeg
            avgData = (avgData * 0.248 * xs)/0.0082 

            avgData /= 1000 

            ax.plot(xs,avgData,color=colorDict[solvent],linewidth=2) 


        cbar_ax = fig2.add_axes([0.80, 0.15, 0.03, 0.7])
        cbar = fig2.colorbar(CS3, cax=cbar_ax,ticks=[0+(binSize/10)/2,1250-(binSize/10)/2])
        cbar.ax.set_yticklabels(["0-%i ns"%(binSize/10),"%i-1250 ns"%(1250-(binSize/10))]) 

        fig2.savefig('figures/cd_spectra_%s_time_resolved_%.1f.png'%(param,bandwidth),format='png') 
        plt.close(fig2) 

        ax.text(0.97,0.95,r"\textsf{%c}"%indexToLetter[combCol+1],va='top',ha='right',transform=ax.transAxes,fontsize=12)

    ax.set_xlim([190,250]) 
    ax.set_ylim([-35,40]) 
    fig.savefig('figures/combined_cd_spectra_%.1f.png'%(bandwidth),format='png')
    plt.close(fig) 





