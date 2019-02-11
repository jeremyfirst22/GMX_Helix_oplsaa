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

for solvent in 'not_bound_sam' ; do #'water' 'tert' 'not_bound_sam' ; do 
#for solvent in 'tert' ; do 
    logFile=${TOP}/$solvent/good_turing/good_turing.log
    errFile=${TOP}/$solvent/good_turing/good_turing.err

    if [ $solvent == 'water' ] ; then 
        equilTime=600 ##ns
        samplingFol=13 
        samplingUnf=6  ##These are taken directly from the individual GT runs. 
    elif [ $solvent == 'tert' ] ; then 
        equilTime=150 
        samplingFol=5
        samplingUnf=7
    elif [ $solvent == 'not_bound_sam' ] ; then 
        equilTime=150 
        samplingFol=13
        samplingUnf=21
    fi 

    echo $solvent 
    cd $solvent

    if [ ! -d good_turing ] ; then mkdir good_turing ; fi 
    cd good_turing 

    if [ ! -f folded.gro ] ; then 
        echo 'Backbone Backbone' | gmx trjconv -s ../folded/Production/folded_${solvent}.tpr \
            -f ../folded/Production/system_npt.gro \
            -pbc mol \
            -ur compact \
            -center \
            -o folded.gro >> $logFile 2>> $errFile 
    fi 
    check folded.gro 

    if [ ! -f folded.xtc ] ; then 
        echo 'Backbone Backbone' | gmx trjconv -s ../folded/Production/folded_${solvent}.tpr \
            -f ../folded/Production/folded_${solvent}.xtc \
            -pbc mol \
            -ur compact \
            -center \
            -b $((equilTime * 1000)) \
            -dt $((samplingFol * 1000)) \
            -o folded.xtc >> $logFile 2>> $errFile 
    fi 
    check folded.xtc 

    if [ ! -f unfolded.xtc ] ; then 
        echo 'Backbone Backbone' | gmx trjconv -s ../unfolded/Production/unfolded_${solvent}.tpr \
            -f ../unfolded/Production/unfolded_${solvent}.xtc \
            -pbc mol \
            -ur compact \
            -center \
            -b $((equilTime * 1000)) \
            -dt $((samplingUnf * 1000)) \
            -o unfolded.xtc >> $logFile 2>> $errFile 
    fi 
    check unfolded.xtc 

    if [ ! -s both.xtc ] ; then 
        gmx trjcat -f folded.xtc unfolded.xtc \
            -cat \
            -o both.xtc >> $logFile 2>> $errFile 
    fi 
    check both.xtc 

    if [ ! -s rmsd.xpm ] ; then
        echo "Backbone Backbone" | gmx rms -s folded.gro \
            -f both.xtc \
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
        


    
