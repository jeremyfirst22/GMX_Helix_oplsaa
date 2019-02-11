#!/bin/bash

TOP=${PWD}

gtScripts=$TOP/good-turing_scripts

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
    logFile=${TOP}/$solvent/good_turing/good_turing.log
    errFile=${TOP}/$solvent/good_turing/good_turing.err

    if [ $solvent == 'water' ] ; then 
        equilTime=600 ##ns
    else : 
        equilTime=150 
    fi 

    echo $solvent 
    cd $solvent

    if [ ! -d good_turing ] ; then mkdir good_turing ; fi 
    cd good_turing 

    sampling=100 
    if [ ! -s both.xtc ] ; then 
        gmx trjcat -f ../folded/nopbc/nopbc.xtc ../unfolded/nopbc/nopbc.xtc \
            -cat \
            -b $((equilTime*1000)) \
            -dt $sampling \
            -o both.xtc >> $logFile 2>> $errFile 
    fi 
    check both.xtc 

    if [ ! -s rmsd.xpm ] ; then
        echo "Backbone Backbone" | gmx rms -s ../folded/nopbc/nopbc.gro \
            -f both.xtc \
            -dt $sampling \
            -b $equilTime \
            -m rmsd.xpm >> $logFile 2>> $errFile
    fi
    check rmsd.xpm

    if [ ! -s rmsd.dat ] ; then 
        python $gtScripts/xpm2dat.py rmsd.xpm 
    fi 
    check rmsd.dat 

    if [ ! -f good_turing.rmsd.tar ] ; then 
        ##https://github.com/pkoukos/GoodTuringMD
        ##  Adapted to avoid interacitve file name input
        echo "source('$gtScripts/Good_Turing.R')" > run_gt_stat.R

        RScript --no-save --no-restore --verbose run_gt_stat.R >> $logFile  2>> $errFile 
    fi 
    check good_turing.rmsd.tar

    cd ../../
done 
        


    
