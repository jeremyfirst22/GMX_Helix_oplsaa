#!/usr/bin/env python 

import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc_file


rcFile = 'rc_files/paper.rc' 
exp_data = 'Exp_data/cd.txt' 

data = np.genfromtxt(exp_data) 

rc_file(rcFile) 

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(left=0.15,right=0.78, bottom=0.12,top=0.95) 
fig.text(0.474,0.04, r"Wavelength (nm)", ha='center', va='center') 
#fig.text(0.03,0.5, r"Num waters within %s A of Q61 s.c. and %s A of O1G"%(dist, dist) , ha='center', va='center',rotation='vertical') 
fig.text(0.03,0.535, r"Ellipticity on SAM $\thetaup$ (mdeg)", ha='center', va='center',rotation='vertical') 
fig.text(0.99,0.535, r"Molar ellipticity in solution $\thetaup_m$ ", ha='right', va='center',rotation=270) 
fig.text(0.90,0.535, r"(mdeg cm$^2$ dmol$^{-1}$)", ha='left', va='center',rotation=270) 


ax2 = ax.twinx() 

ax2.plot(data[:,0],data[:,1],'b',linewidth=2) 
ax2.plot(data[:,0],data[:,2],'k',linewidth=2) 
ax.plot(data[:,0],data[:,3],'g',linewidth=2) 

##Baseline correct SAM spectra for BeStSel deconvolution
data[:,3] -= data[0,3]  ##First element is abs @ 250 nm. Set this to zero. 
##RedShift correction (blue-shift by 2 nm) 
data[:,3] = np.insert(data[:-1,3],0,0)

tbut = np.array([data[:,0],data[:,1] ])  ##20 residues. Molar ellipticity -> mean residue ellipticity
water = np.array([data[:,0],data[:,2]  ])  ##20 residues. Molar ellipticity -> mean residue ellipticity
sam = np.array([data[:,0],data[:,3]  ] )  ##ellipticity

np.savetxt('Exp_data/tbut.txt',tbut.T,fmt='%3.6f') 
np.savetxt('Exp_data/water.txt',water.T,fmt='%3.6f') 
np.savetxt('Exp_data/sam.txt',sam.T,fmt='%3.6f') 


ax2.set_ylim([-18,0]) 
ax.set_ylim([-1,2.5]) 
ax.set_xlim([200,250]) 
#ax2.set_xlabel("Molar Ellipticity on SAM [$\Theta$] (mdeg)") 



plt.savefig('figures/nachos_cd_spectra.png',format='png') 
