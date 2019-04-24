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

for solvent in 'water' 'tert' 'sam' 'single_bound_sam' 'not_bound_sam' 'disordered_sam' 'disordered_single_bound_sam' 'disordered_not_bound_sam' ; do 
#for solvent in 'tert' ; do 
    logFile=${TOP}/$solvent/sham/sham.log
    errFile=${TOP}/$solvent/sham/sham.err

    if [ $solvent == 'water' ] ; then 
        equilTime=600 ##ns
    else : 
        equilTime=150 
    fi 
    #equilTime=0

    echo $solvent 
    cd $solvent

    if [ ! -d sham ] ; then mkdir sham ; fi 
    cd sham 

    cat ../folded/rgyr/gyrate.xvg > data.xvg 
    cat ../unfolded/rgyr/gyrate.xvg >> data.xvg 
    echo '&' >> data.xvg 
    cat ../folded/rmsd/rmsd.xvg >> data.xvg 
    cat ../unfolded/rmsd/rmsd.xvg >> data.xvg 

    gmx sham -f data.xvg \
        -ls gibbs.xpm \
        -lp prob.xpm \
        -n 2 \
        -b $((equilTime*1000)) \
        -xmin 0.75 0.150 0 \
        -xmax 1.75 1.100 0 \
        -nlevels 25  \
        -ngrid 100 >> $logFile 2>> $errFile 
    check gibbs.xpm 

    if [ -f gibbs.dat ] ; then 
        rm gibbs.dat 
    fi 

    python ../../analysis_scripts/xpm2dat.py gibbs.xpm 
    check gibbs.dat 

    cd ../../
done 
        


    
