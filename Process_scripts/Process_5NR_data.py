import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc_file

sampleList = ["5", "15"] 
for sample in sampleList : 
    ###
    #  Figure set up and labeling
    ###
    rc_file('rc_files/paper.rc')
    
    left, right = 0.15, 0.95
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
    #fileList = ["NR_data/ORNL15_air_R.txt", 
    envList = ["air", "D2O", "7030"] 
    #envList = ["D2O"] 
    env2Color = {"air":'black', "D2O":'green', "7030":'blue'}
    env2Stag = {"air":1.5, "D2O":0, "7030":-1.5}
    env2label= {"air":"Air", "D2O":r"D$_2$O", "7030":"70/30:H$_2$O/D$_2$O"}
    
    for env in envList : 
        rFileName = "NR_data/ORNL%s_%s_R.txt"%(sample,env)
        mFileName = "NR_data/fit_ORNL%s_%s.txt"%(sample,env)
        print rFileName
    
        try : 
            rData = np.genfromtxt(rFileName) 
        #mData = np.genfromtxt(mFileName) 
        except IOError : 
            print "Warning: %s file not found. Skipping"%env 
            break
    
        qValues = rData[:,0] 
        #refl = np.log10(rData[:,1]) 
        refl = rData[:,1]
        refl *= 10**env2Stag[env]
    
        reflErr = rData[:,2]  
        reflErr *= 10**env2Stag[env]
        
        ax1.scatter(qValues, refl,s=5.0,marker='o', color='none', edgecolor=env2Color[env],linewidth=0.5) 
        ax1.errorbar(qValues,refl, yerr=reflErr,capsize=0,color=env2Color[env],elinewidth=0.5,linewidth=0 ) 
    
        ax1.set_yscale('log') 
    
        try : 
            mData = np.genfromtxt(mFileName,skip_header=1)
        #mData = np.genfromtxt(mFileName) 
        except IOError : 
            print "Warning: %s file not found. Skipping"%env 
            continue 
        mData[:,1] += env2Stag[env]
        ax1.plot(mData[:,0],10**(mData[:,1]),color=env2Color[env]) 
    
    ax1.legend(loc=1)
    
    ###
    #  Subplot 1: SLD Plots
    ### 
    sldData = np.genfromtxt('NR_data/sld_graph_ORNL%s.txt'%sample,skip_header=1) 
    
    ax2.plot(sldData[:,0]-07, sldData[:,2], color=env2Color["air"]) 
    ax2.plot(sldData[:,0]-00, sldData[:,3], color=env2Color["D2O"]) 
    ax2.plot(sldData[:,0]-11, sldData[:,4], color=env2Color["7030"]) 
    
    airLayer = 480 
    ax2.axvline(airLayer, color='k', linestyle='--')
    ax2.axvline(airLayer - 20.8, color='k', linestyle='--') 
    ax2.axvline(airLayer - 20.8 - 10.2, color='k', linestyle='--') 
    
    #ax2.axvline()
    #ax2.axvline()
    
    ax2.set_xlim([400,488]) 
    ax2.set_ylim([0,6.2]) 
    
    
    
    
    
    
    
    fig.savefig('figures/NR_data_ORNL%s.png'%sample,format='png')
