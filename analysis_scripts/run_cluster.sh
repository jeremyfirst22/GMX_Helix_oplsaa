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
    logFile=${TOP}/$solvent/cluster/cluster.log
    errFile=${TOP}/$solvent/cluster/cluster.err

    if [ $solvent == 'water' ] ; then 
        equilTime=600 ##ns
    else : 
        equilTime=150 
    fi 

    echo $solvent 
    cd $solvent

    if [ ! -d cluster ] ; then mkdir cluster ; fi 
    cd cluster 

    if [ ! -s both.xtc ] ; then 
        gmx trjcat -f ../folded/nopbc/nopbc.xtc ../unfolded/nopbc/nopbc.xtc \
            -cat \
            -b $((equilTime*1000)) \
            -dt 10000 \
            -o both.xtc >> $logFile 2>> $errFile 
    fi 
    check both.xtc 

    if [ ! -f clusters.pdb ] ; then 
        echo "Backbone System" | gmx cluster -f both.xtc \
            -s ../folded/nopbc/nopbc.gro \
            -dt 100 \
            -cutoff 0.25 \
            -method gromos \
            -tr clust-trans.xpm \
            -cl clusters.pdb \
            -sz clust-size.xvg >> $logFile 2>> $errFile 
    fi 
    check clust-size.xvg clusters.pdb 

    cd ../../
done 
        


    
