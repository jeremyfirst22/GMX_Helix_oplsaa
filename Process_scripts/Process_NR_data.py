import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc_file


envList = ["Air", "D2O", "7030"] 
sampleNames = {
        "2" : {"D2O": 390495, "Air": 229450 ,"7030": 390496},  
        "5" : {"D2O": 229453, "Air": 229449, "7030": 229452}, 
        "15" : {"D2O": 390516, "Air": 229451,"7030": 390520}  
        }
chi2_values = {
        "2": {"Air":1.8, "D2O":2.0, "7030":1.9}, 
        "5": {"Air":2.3, "D2O":1.6, "7030":4.0}, 
        "15":{"Air":1.8, "D2O":1.3, "7030":1.9} 
        }

#env2Chi_2 = {"Air":1.8, "D2O":2.0, "7030":2.0}
#env2Chi_5  = {"Air":2.4, "D2O":1.3, "7030":1.9}
#env2Chi_15 = {"Air":1.9, "D2O":1.2, "7030":1.8}

env2Color = {"Air":'black', "D2O":'green', "7030":'blue'}
env2Stag = {"Air":1.5, "D2O":0, "7030":-1.5}
env2label= {"Air":"Air", "D2O":r"D$_2$O", "7030":"70/30 D$_2$O/H$_2$O"}

sampleList = ["2", "5" , "15"] 
for sample in sampleList : 
    ###
    #  Figure set up and labeling
    ###
    rc_file('rc_files/paper.rc')
    
    left, right = 0.13, 0.95
    bottom, top = 0.10, 0.98
    hspace, wspace = 0.25, 0.00
    
    fig,axarr = plt.subplots(2,1,figsize=(3.25,4))
    fig.subplots_adjust(left=left, right=right, top=top,bottom=bottom,hspace=hspace, wspace=wspace)
    fig.text(0.01,top-(top-bottom-hspace/2)/4,"Log(reflectivity)", rotation='vertical', ha='left',va='center')
    fig.text(0.01,bottom+(top-bottom-hspace/2)/4,"SLD ($\\times 10^{-6}$ $\AA^{-2}$)", rotation='vertical', ha='left',va='center')
    
    fig.text((right-left)/2+left, 0.50, "$Q$ ($\AA^{-1}$)", ha='center', va='bottom')
    fig.text((right-left)/2+left, 0.00, "Distance from Si interface ($\AA$)", ha='center', va='bottom')
    
    ax1, ax2 = axarr
    
    ax1.text(0.05,0.95,r"\textsf{A}",va='top',ha='left',transform=ax1.transAxes,fontsize=12)
    ax2.text(0.05,0.95,r"\textsf{B}",va='top',ha='left',transform=ax2.transAxes,fontsize=12)
    
    ###
    #  Subplot 1: Reflectivity profile
    ### 
    for env in envList : 
        rFileName = "NR_data_2/reflectivity_fits/%s/__model-refl.dat"%(sampleNames[sample][env]) 
#        mFileName = "NR_data/reflectivitity_data/%s/__model-refl.dat"%(sample,env)
        print rFileName
    
        try : 
            rData = np.genfromtxt(rFileName, skip_header=3) 
        #mData = np.genfromtxt(mFileName) 
        except IOError : 
            print "Warning: %s file not found. Skipping"%env 
            continue 
    
        qValues = rData[:,0] 
        #refl = np.log10(rData[:,1]) 
        refl = rData[:,2]
        refl *= 10**env2Stag[env]
    
        reflErr = rData[:,3]  
        reflErr *= 10**env2Stag[env]

        mData = rData[:,4] 
        mData *= 10**env2Stag[env]
        
        ax1.scatter(qValues, refl,s=5.0,marker='o', color='none', edgecolor=env2Color[env],linewidth=0.5,zorder=10) 
        ax1.errorbar(qValues,refl, yerr=reflErr,capsize=0,color=env2Color[env],elinewidth=0.5,linewidth=0,zorder=5 ) 
    
        ax1.set_yscale('log') 

        ax1.text(0.240, 10**-6 * 10**env2Stag[env], "$\chi^2 = %.1f$"%chi2_values[sample][env], color=env2Color[env], ha='right', va='center',fontsize='small') 
    
        ax1.plot(qValues,mData,color=env2Color[env],zorder=1,linewidth=0.5) 
    
#    ax1.legend(loc=1)#,fontsize=6)
    ax1.set_ylim([10**-9, 10**4]) 
    ax1.set_xlim(-0.005, 0.245) 
#    ax1.set_yticks(np.geomspace(10**-9, 10**4, num=16))  
    ax1.set_yticklabels(np.log10(ax1.get_yticks()) ) 
    print ax1.get_yticks() 
    
    ###
    #  Subplot 2: SLD Plots
    ### 
    sample2offset = {"2":400, "15":480, "5":200}
    sample2window= {"2":[400,480], "15":[400,488], "5":[210,310]}
    
    alignLayer = 1    ## SiOx layer
    ###
    #  Read through all environments and find avg interface and layer thicknessess
    ###
    avgEnv, avgPepT, avgSamT, index = 0,0,0,0
    for env in envList : 
        slabFileName = "NR_data_2/reflectivity_fits/%s/__model-slabs.dat"%(sampleNames[sample][env]) 
        try : 
            slabData = np.genfromtxt(slabFileName, skip_header=1) 
        except IOError : 
            print "Skipping SLD"
            continue 
        index += 1
        thickValues = slabData[:,0]
        ###
        #  Apparently the x-axis of all the Air data is inverted for no reason at all. 
        ###
        if not env == "Air" : 
            thickValues = np.flip(thickValues) 

#        zSilicon = np.sum(thickValues[:-1]) 
        envI = 1 ##First layer (counting from top) 
        avgEnv  += np.sum(thickValues[:-envI])           ## Position of Env interface w.r.t. bottom of SiOx
        avgPepT += thickValues[5]                     ## Thickness of peptide layer
        avgSamT += thickValues[4]                     ## Thickness of SAM layer
    avgEnv /= index
    avgPepT /= index
    avgSamT /= index
    
    print "Peptide layer height = ", avgPepT, "SAM layer height = ", avgSamT

    ###
    #  Read through again and plot each profile, aligned by Au/SAM interface
    ### 
    for env in envList : 
        sldFileName = "NR_data_2/reflectivity_fits/%s/__model-profile.dat"%(sampleNames[sample][env]) 
        print sldFileName
        try : 
            sldData = np.genfromtxt(sldFileName, skip_header=1) 
        except IOError : 
            print "Skipping SLD"
            continue 

        slabFileName = "NR_data_2/reflectivity_fits/%s/__model-slabs.dat"%(sampleNames[sample][env]) 
        try : 
            slabData = np.genfromtxt(slabFileName, skip_header=1) 
        except IOError : 
            print "Skipping SLD"
            continue 

        zValues   = sldData[:,0]
        sldValues = sldData[:,1]
        thickValues = slabData[:,0] 

        ###
        #  Apparently the x-axis of liquid environments is opposite that off air (Si layer at max height instead of set to zero)
        #       1. Translate so that Silicon interface is set to zero
        #       2. Flip axis 
        ###
        if not env == "Air" : 
            zValues -= np.sum(thickValues[:]) 
            zValues *= -1

        ###
        #  Align by average d of Au interface
        ###
        zEnv = np.sum(thickValues[:-envI])
        zValues -= zEnv - avgEnv
        zEnv -= zEnv - avgEnv
    
        ###
        #  Plot profiles 
        ### 
        ax2.plot(zValues, sldValues , color=env2Color[env],label=env2label[env]) 
        #ax2.axvline(zEnv, color='gray', linestyle='--')

    ax2.axvline(avgEnv, color='k', linestyle='--')
    ax2.axvline(avgEnv - avgPepT, color='k', linestyle='--')
    ax2.axvline(avgEnv - avgPepT - avgSamT, color='k', linestyle='--')

    if sample in ["15"] : rightCushion = 10
    else :                rightCushion = 20   ##Profile plots cut early for Sample 15

    leftCushion = 45
    xMax = avgEnv + rightCushion 
    xMin = avgEnv - avgPepT - avgSamT - leftCushion 
    ax2.set_xlim([xMin,xMax]) 
#    ax2.set_ylim([0,6.2]) 

    ###
    #  Label regions of interest; Au: Gold layer; 1: SAM layer; 2: Peptide layer; 3: Environment
    ###
    ax2.text(xMax - (xMax - avgEnv)/2,   1.5, "3", ha='center', va='center') 
    ax2.text(avgEnv - avgPepT/2,   1.5, "2", ha='center', va='center') 
    ax2.text(avgEnv - avgPepT - avgSamT/ 2,   1.5, "1", ha='center', va='center') 
    ax2.text(xMin + (avgEnv - avgPepT - avgSamT - xMin)/2,   1.5, "Au", ha='center', va='center') 

    
    fig.legend(loc='upper right', bbox_to_anchor=(right,top)) 

    
    fig.savefig('figures/NR_data_%s.png'%sample,format='png')

