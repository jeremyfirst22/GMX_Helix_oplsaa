#!/bin/bash

TOP=${PWD}

check(){
   for arg in $@ ; do
        if [ ! -s $arg ] ; then
            echo ; echo "ERROR: $arg missing. Exitting"
            exit
            fi
        done
}

for solvent in 'water' 'tert' 'not_bound_sam' ; do 
#for solvent in 'tert' ; do 
    logFile=${TOP}/$solvent/rgyr_rmsd_cluster/log.log
    errFile=${TOP}/$solvent/rgyr_rmsd_cluster/log.err

    echo $solvent 
    cd $solvent

    if [ ! -d rgyr_rmsd_cluster ] ; then mkdir rgyr_rmsd_cluster ; fi 
    cd rgyr_rmsd_cluster 

    if [ ! -f rmsd.xvg ] ; then 
        echo "Backbone Backbone" | gmx rms -f ../cluster/clusters.pdb \
            -s ../../StartingStructures/folded.pdb \
            -fit rot+trans \
            -o rmsd.xvg >> $logFile 2>> $errFile 
    fi 
    check rmsd.xvg 

    if [ ! -f gyrate.xvg ] ; then 
        echo "Protein" | gmx gyrate -f ../cluster/clusters.pdb \
            -s ../cluster/clusters.pdb \
            -o gyrate.xvg >> $logFile 2>> $errFile 
    fi 
    check gyrate.xvg 

    cd ../../
done 
        


    
