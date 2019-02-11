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

for solvent in 'water' 'tert' 'sam' ; do 
#for solvent in 'tert' ; do 
    logFile=${TOP}/$solvent/dssp_unfolding/dssp_unfolding.log
    errFile=${TOP}/$solvent/dssp_unfolding/dssp_unfolding.err

    echo $solvent 
    cd $solvent

    if [ ! -d dssp_unfolding ] ; then mkdir dssp_unfolding ; fi 
    cd dssp_unfolding 

    if [ ! -s unfolding.xtc ] ; then 
        gmx trjcat -f ../prep/Heating/heating.xtc ../prep/Heated_nvt/heated_nvt.xtc ../prep/Cooling/cooling.xtc \
            -cat \
            -o unfolding.xtc >> $logFile 2>> $errFile 
    fi 
    check unfolding.xtc 

    if [ ! -f scount.xvg ] ; then
        echo 'Protein' | gmx do_dssp -f unfolding.xtc \
            -s ../prep/Heating/heating.tpr \
            -ver 1 \
            -sss HGI \
            -ssdump ssdump.dat \
            -o ss.xpm \
            -a area.xpm \
            -ta totarea.xvg \
            -aa averarea.xvg \
            -sc scount.xvg >> $logFile 2>> $errFile
        fi
    check scount.xvg ss.xpm area.xpm

    if [ ! -f ss.png ] ; then 
        gmx xpm2ps -f ss.xpm \
            -di $TOP/m2p_files/unfolding.m2p \
            -by 10 \
            -o ss.eps >> $logFile 2>> $errFile 

        ps2pdf ss.eps ss.pdf 
        convert ss.pdf ss.png 
    fi 
    check ss.png

    cd ../../
done 
        


    
