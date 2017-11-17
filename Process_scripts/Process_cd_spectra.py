#!/usr/bin/env python 

projectDir='Users/jeremyfirst/GMX_Helix_oplsaa'
#rcFile='paper.rc'
rcFile='presentation.rc'
saveDir='figures'

import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc_file
import glob as glob 
import os 
import sys 

rc_file('rc_files/%s'%rcFile) 

binSize= 500

top=os.getcwd()

colorDict = {'water':'k','tert':'b','sam':'g'}

for param in [ 'hirst','woody' ] : 
    totFig = plt.figure() 
    totAx = totFig.add_subplot(111) 
    for solvent in ['water', 'tert', 'sam'] : 
        totAvg = 0 
        color = colorDict[solvent] 
        for state in ['folded', 'unfolded'] : 
            fig = plt.figure() 
            ax = fig.add_subplot(111) 
            plt.ylabel(r"[$\theta$] (10$^3$ deg cm$^2$ dmol$^{-1}$)",ha='center',va='center',rotation='vertical')  
            plt.xlabel(r"$\lambda$ (nm)",ha='center',va='center') 
            plt.xlim([190,250]) 

            datafiles = glob.glob("%s/%s/cd_spectra/output-%s-%s-%s/frame*.cd"%(solvent,state,solvent,state,param) ) 
            if len(datafiles) == 0 :  
                print "No spectra found for %s %s with %s parameters" %(solvent,state,param)
                continue 

            done = False 
            index = 0 
            try : 
                data = np.genfromtxt("%s/%s/cd_spectra/output-%s-%s-%s/frame%i.cd"%(solvent,state,solvent,state,param,index)) 
            except : 
                print "Processed %i for %s %s %s"%(index, solvent, state, param)
                continue 
            avgdata = np.zeros_like(data[:,1]) 
            datax = data[:,0]
            
            for j in np.arange(0,len(datafiles)/binSize) : 
                bindata = np.zeros_like(data[:,1]) 
                index = 0 
                for i in np.arange(j*binSize,(j+1)*binSize) : 
                    try : 
                        data = np.genfromtxt("%s/%s/cd_spectra/output-%s-%s-%s/frame%i.cd"%(solvent,state,solvent,state,param,i)) 
                    except : 
                        print "Processed %i for %s %s %s"%(i, solvent, state, param)
                        break
                    bindata += data[:,1]
                    avgdata += data[:,1]
                    index +=1
                bindata /= index * 1000 
                label="%i-%i ns"%(binSize*j/10,binSize*(j+1)/10) 
                plt.plot(datax,bindata,label=label) 
                continue 
            plt.title("From %s in %s"%(state,solvent)) 
            plt.ylim([-25,25]) 

            colormap = plt.cm.jet
            colors = [colormap(i) for i in np.linspace(0,1,len(ax.lines))]
            for i,j in enumerate(ax.lines) : 
                j.set_color(colors[i]) 

            plt.legend() 
            plt.plot(datax,np.zeros_like(datax),linestyle='--',color='k') 
            fig.savefig("figures/cd_spectra_%s_%s_%s.pdf"%(state,solvent,param),format='pdf' ) 
            print "Figure complete for %s %s using %s"%(state,solvent, param) 
            fig.clf()   

            fig = plt.figure() 
            ax = fig.add_subplot(111) 
            ax.set_ylabel(r"[$\theta$] (10$^3$ deg cm$^2$ dmol$^{-1}$)",ha='center',va='center',rotation='vertical')  
            ax.set_xlabel(r"$\lambda$ (nm)",ha='center',va='center') 
            ax.set_title("From %s in %s"%(state,solvent)) 
            ax.set_xlim([190,250]) 
            ax.plot(datax,avgdata/(len(datafiles) * 1000) ,linewidth=5,label='Trajectory average',color=color ) 
            ax.set_ylim([-25,25]) 
            ax.plot(datax,np.zeros_like(datax),linestyle='--',color='k') 
            #plt.legend() 
            fig.savefig("figures/cd_spectra_avg_%s_%s_%s.pdf"%(state,solvent,param),format='pdf' ) 
            print "Avg. figure complete for %s %s using %s"%(state,solvent, param) 
            fig.clf()   
            plt.close() 

            try : 
                totAvg += avgdata
            except : 
                totAvg = np.zeros_like(avgdata) 
                totAvg += avgdata

        print len(datax), len(totAvg) 
        totAx.plot(datax,totAvg/(2* len(datafiles) * 1000) ,linewidth=3,label=solvent,color=color ) 
    totAx.plot(datax,np.zeros_like(datax),linestyle='--',color='k') 
    totAx.set_title("Avg calculated CD from all states")
    totAx.set_xlim([190,250]) 
    totAx.set_ylim([-10,15]) 
    totAx.legend() 
    totFig.savefig('figures/combined_cd_spectra_%s.pdf'%param,format='pdf')  


                 


