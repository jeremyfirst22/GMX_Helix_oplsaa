#!/bin/bash

#dim=5.964
numMols=12
spacing=0.497 ##nm 
glyDist=0.6822 ##nm from geometry optimization, Spartan '16, DFT wB97X-D, 6-311+g**
glyRest=1000   #kJ/mol/nm
sulRest=200000 #kJ/mol/nm

usage(){
    echo "USAGE: $0 <PDB file {molec.pdb} > " 
    exit
} 

if [ -z $1 ] ; then 
    usage 
    fi 

fileName=$1
if [ ! -f $fileName ] ; then 
    echo "ERROR: $fileName not found " 
    exit 
    fi 
if [[ $fileName == *.pdb ]] ; then 
    MOLEC=$(basename $fileName) 
    MOLEC=${MOLEC%.*}_sam
else 
    echo "ERROR: Input file must be PDB file (*.pdb)" 
    exit 
    fi 
if [ ! -d $MOLEC ] ; then mkdir $MOLEC ; fi 
if [ ! -f $MOLEC/$fileName ] ; then cp $fileName $MOLEC/$MOLEC.pdb ; fi 

TOP=${PWD}
MDP=$TOP/mdp_files
logFile=$TOP/$MOLEC/$MOLEC.log
errFile=$TOP/$MOLEC/$MOLEC.err
FF=$TOP/GMXFF
forceField=oplsaa
if [ ! -d $FF/$forceField.ff ] ; then 
    echo ; echo "ERROR: FF not found" 
    exit
    fi 

if [ -f StartingStructures/decanethiol.pdb ] ; then 
    cp StartingStructures/decanethiol.pdb $MOLEC/. 
else 
    echo ; echo "ERROR: StartingStructures/decanethiol.pdb" 
    exit
    fi 

check(){
    for arg in $@ ; do 
        if [ ! -s $arg ] ; then
            echo "ERROR: $arg missing. Exitting" 
            exit 
            fi 
        done 
} 

clean(){
    if [ -d $forceField.ff ] ; then rm -r $forceField.ff *.dat ; fi 
}

create_dir(){
    if [ -z $1 ] ; then 
        echo "ERROR: create_dir requires argument. " ; exit ; fi  

    dirName=$1 
    if [ ! -d $dirName ] ; then mkdir $dirName ; fi  
    
    if [ ! -d $dirName/$forceField.ff ] ; then 
           if [ -d $FF/$forceField.ff ] ; then 
               cp -r $FF/$forceField.ff $dirName
               cp $FF/*.dat $dirName 
           else 
               echo "FF not found" 
               exit 
               fi  
           fi  
}

build_SAM(){
    printf "\t\tBuilding SAM layer........................" 
    if [ ! -f Build_SAM/bottom.gro ] ; then 
        create_dir Build_SAM 
        
        cp $MOLEC.pdb Build_SAM/. 
        cp decanethiol.pdb Build_SAM/. 
        cd Build_SAM

        echo 'LIG' | gmx editconf -f decanethiol.pdb \
            -princ \
            -o oriented.gro >> $logFile 2>> $errFile 
        check oriented.gro 

        gmx editconf -f oriented.gro \
            -rotate 180 90 0 \
            -o upright.gro >> $logFile 2>> $errFile 
        check upright.gro 

        gmx editconf -f upright.gro \
            -rotate 30 0 0 \
            -center 0 0 0 \
            -o inplane.gro >> $logFile 2>> $errFile  
        check inplane.gro 

        gmx editconf -f inplane.gro \
            -rotate 0 22 0 \
            -center 0 0 0 \
            -o twisted.gro >> $logFile 2>> $errFile 
        check twisted.gro 

        echo "#!/usr/bin/env python
import math as m

numMols=$numMols
spacing=$spacing
step=spacing * m.cos(m.pi / 6 ) 

for i in range(numMols) : 
    for j in range(numMols) : 
        x,y,z = i*step, j*spacing, 0.0 
        if not i%2 == 0 : 
             y += spacing / 2
        print \"%.3f  %.3f  %.3f\"%(x,y,z) 
" > make_position.py 

        python make_position.py > position.dat 
        check position.dat 
        
        xdim=`echo "$numMols * $spacing * c(4*a(1) / 6) " | bc -l`
        ydim=`echo "$numMols * $spacing " | bc -l`
        zdim=$xdim

        gmx insert-molecules -ci twisted.gro \
            -ip position.dat \
            -rot none \
            -box $xdim $ydim $zdim \
            -o layer.gro >> $logFile 2>> $errFile 
        check layer.gro 

        gmx editconf -f layer.gro \
            -bt triclinic \
            -o boxed.gro >> $logFile 2>> $errFile  
        check boxed.gro 

        gmx pdb2gmx -f boxed.gro \
            -p layer.top \
            -ff $forceField \
            -water tip3p \
            -o boxed.gro >> $logFile 2>> $errFile 
        check layer.top 

        includeLine=`cat -n layer.top | grep '; Include Position restraint file' | awk '{print $1}' | tail -n1 `
        head -n $includeLine layer.top > decanethiol.top 
        echo '#ifdef POSSULRES 
#include "posre_SUL.itp" 
#endif' >> decanethiol.top 
        ((includeLine++)) 
        tail -n +$includeLine layer.top >> decanethiol.top

        echo "[ position_restraints ]" > posre_SUL.itp 
        echo ";; Pin sulfur atoms to spacing found on SAM" >> posre_SUL.itp 
        for atom in `grep LIG boxed.gro | grep S1 | awk '{print $3}'` ; do 
            printf "%6i%6i%10.f%10.f%10.f\n" $atom 1 $sulRest $sulRest $sulRest >> posre_SUL.itp 
        done 

        zshift=`grep LIG boxed.gro | grep " H22 " | awk '{print $6}' | sort -n | uniq | head -n1`
        yshift=`grep LIG boxed.gro | grep " H22 " | awk '{print $5}' | sort -n | uniq | head -n1`
        xshift=`grep LIG boxed.gro | grep " H22 " | awk '{print $4}' | sort -n | uniq | head -n1`
    
        zshift=`echo "-$zshift + 0.10" | bc -l | awk '{printf "%f", $0}'`
        yshift=`echo "-$yshift " | bc -l | awk '{printf "%f", $0}'`
        xshift=`echo "-$xshift " | bc -l | awk '{printf "%f", $0}'`

        gmx editconf -f boxed.gro \
            -translate $xshift $yshift $zshift \
            -o bottom.gro >> $logFile 2>> $errFile 
        check bottom.gro

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
    check Build_SAM/posre_SUL.itp Build_SAM/decanethiol.top Build_SAM/bottom.gro 
} 

layer_relax(){
    printf "\t\tSteep and NVT relax of SAM layer.........." 
    if [ ! -f Relax_SAM/nvt_relax.nopbc.gro ] ; then 
        create_dir Relax_SAM
        
        cp Build_SAM/decanethiol.top Relax_SAM/.
        cp Build_SAM/bottom.gro Relax_SAM/.
        cp Build_SAM/*.itp Relax_SAM/. 
        cd Relax_SAM
        
        gmx grompp -f $MDP/sam_steep.mdp \
            -c bottom.gro \
            -p decanethiol.top \
            -o steep.tpr >> $logFile 2>> $errFile 
        check steep.tpr 

        gmx mdrun -deffnm steep >> $logFile 2>> $errFile 
        check steep.gro 

        gmx grompp -f $MDP/sam_nvt_relax.mdp \
            -c steep.gro \
            -p decanethiol.top \
            -o nvt_relax.tpr >> $logFile 2>> $errFile 
        check nvt_relax.tpr 

        if [ ! -f nvt_relax.gro ] ; then 
            gmx mdrun -deffnm nvt_relax >> $logFile 2>> $errFile 
            fi 
        check nvt_relax.gro 

        echo '0' | gmx trjconv -f nvt_relax.gro \
            -s nvt_relax.tpr \
            -pbc nojump \
            -o nvt_relax.nopbc.gro >> $logFile 2>> $errFile 
        check nvt_relax.nopbc.gro

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

protein_steep(){
    printf "\t\tProtein steep............................." 
    if [ ! -f Protein_steep/protein_steep.gro ] ; then 
        create_dir Protein_steep
        
        cp $MOLEC.pdb Protein_steep/.
        cd Protein_steep

        ## 3, 3 -- None, None for termini options
        echo '3 3' | gmx pdb2gmx -f $MOLEC.pdb \
            -p $MOLEC.top \
            -ff $forceField \
            -ter \
            -water tip3p \
            -o $MOLEC.gro >> $logFile 2>> $errFile 
        check $MOLEC.gro 

        xdim=`tail -n1 ../Relax_SAM/nvt_relax.gro | awk '{print $1}'`
        ydim=`tail -n1 ../Relax_SAM/nvt_relax.gro | awk '{print $2}'`
        zdim=`tail -n1 ../Relax_SAM/nvt_relax.gro | awk '{print $3}'`

        echo 'Backbone' | gmx editconf -f $MOLEC.gro \
            -box $xdim $ydim $zdim \
            -bt tric \
            -princ \
            -center 0 0 0 \
            -o aligned.gro >> $logFile 2>> $errFile
        check aligned.gro 

        ## 90 deg around z axis puts long axis in long dimension of box
        ## 60 or 30 deg around x axis puts CA of Gly at bottom, with Leus closest to SAM
        if [[ "$MOLEC" == "folded_"* ]] ; then 
            gmx editconf -f aligned.gro \
                -rotate -60 0 90 \
                -o rotated.gro >> $logFile 2>> $errFile 
        elif [[ "$MOLEC" == "unfolded_"* ]] ; then 
            gmx editconf -f aligned.gro \
                -rotate 30 0 90 \
                -o rotated.gro >> $logFile 2>> $errFile 
        else 
            echo "ERROR: I don't know how to orient the protein!" 
            echo "$MOLEC" 
            exit 
            fi 
        check rotated.gro 

        ##This is the distance to the edge of the box. 
        xshift=`echo "$xdim / 2" | bc -l`
        yshift=`echo "$ydim / 2" | bc -l`

        zshift=`grep GLY rotated.gro | grep CA | awk '{print $6}' | sort -nr | tail -n1`
        zshift=`echo "$zshift * -1 + $glyDist" | bc -l`
        gmx editconf -f rotated.gro \
            -translate $xshift $yshift $zshift \
            -o translated.gro >> $logFile 2>> $errFile  
        check translated.gro 

        gmx grompp -f $MDP/protein_steep.mdp \
            -c translated.gro \
            -p $MOLEC.top \
            -o protein_steep.tpr >> $logFile 2>> $errFile 
        check protein_steep.tpr 

        gmx mdrun -deffnm protein_steep >> $logFile 2>> $errFile 
        check protein_steep.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi
} 

solvate(){
    printf "\t\tSolvating protein........................." 
    if [ ! -f Solvate/neutral.top ] ; then 
        create_dir Solvate
        
        cp Protein_steep/protein_steep.gro Solvate/. 
        cp Protein_steep/$MOLEC.top Solvate/. 
        cd Solvate

        xdim=`tail -n1 protein_steep.gro | awk '{print $1}'`
        ydim=`tail -n1 protein_steep.gro | awk '{print $2}'`
        zdim=`tail -n1 protein_steep.gro | awk '{print $3}'`

        gmx solvate -cp protein_steep.gro \
            -box $xdim $ydim $zdim \
            -p $MOLEC.top \
            -o solvated.gro >> $logFile 2>> $errFile 
        check solvated.gro

        gmx grompp -f $MDP/vac_md.mdp \
            -p $MOLEC.top \
            -c solvated.gro \
            -o genion.tpr >> $logFile 2>> $errFile 
        check genion.tpr
        
        echo 'SOL' | gmx genion -s genion.tpr \
            -neutral \
            -nname 'CL' \
            -pname 'NA' \
            -o neutral.gro >> $logFile 2>> $errFile 
        check neutral.gro 

        ## 3, 3 -- None, None for terimini options
        echo '3 3' | gmx pdb2gmx -f neutral.gro \
            -ff $forceField \
            -water tip3p \
            -p neutral.top \
            -ter \
            -o neutral.gro >> $logFile 2>> $errFile 
        check neutral.top 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

solvent_steep(){
    printf "\t\tSolvent steep............................." 
    if [ ! -f Solvent_steep/solvent_steep.gro ] ; then 
        create_dir Solvent_steep
        
        cp Solvate/neutral.gro Solvent_steep/. 
        cp Solvate/neutral.top Solvent_steep/. 
        cp Solvate/*.itp Solvent_steep/. 
        cd Solvent_steep

        gmx grompp -f $MDP/solvent_steep.mdp \
            -p neutral.top \
            -c neutral.gro \
            -o solvent_steep.tpr >> $logFile 2>> $errFile 
        check solvent_steep.tpr 

        gmx mdrun -deffnm solvent_steep >> $logFile 2>> $errFile 
        check solvent_steep.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

solvent_nvt(){
    printf "\t\tSolvent NVT relaxation...................." 
    if [ ! -f Solvent_nvt/solvent_nvt.gro ] ; then 
        create_dir Solvent_nvt
        
        cp Solvent_steep/solvent_steep.gro Solvent_nvt/. 
        cp Solvent_steep/neutral.top Solvent_nvt/. 
        cp Solvent_steep/*.itp Solvent_nvt/. 
        cd Solvent_nvt

        gmx grompp -f $MDP/solvent_nvt_relax.mdp \
            -c solvent_steep.gro \
            -p neutral.top \
            -o solvent_nvt.tpr >> $logFile 2>> $errFile 
        check solvent_nvt.tpr 

        gmx mdrun -deffnm solvent_nvt >> $logFile 2>> $errFile 
        check solvent_nvt.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

solvent_npt(){
    printf "\t\tSolvent NPT semiisotropic relaxation......" 
    if [ ! -f Solvent_npt/solvent_npt.gro ] ; then 
        create_dir Solvent_npt
        
        cp Solvent_nvt/solvent_nvt.gro Solvent_npt/. 
        cp Solvent_nvt/neutral.top Solvent_npt/. 
        cp Solvent_nvt/*.itp Solvent_npt/. 
        cd Solvent_npt

        gmx grompp -f $MDP/solvent_semiisotropic_npt_relax.mdp \
            -c solvent_nvt.gro \
            -p neutral.top \
            -o solvent_npt.tpr >> $logFile 2>> $errFile 
        check solvent_npt.tpr 

        gmx mdrun -deffnm solvent_npt >> $logFile 2>> $errFile 
        check solvent_npt.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

build_system(){
    printf "\t\tBuilding surface+protein+water system....." 
    if [ ! -f Build_system/system.top ] ; then 
        create_dir Build_system 
        
        cp Solvent_npt/solvent_npt.gro Build_system/. 
        #cp Solvent_npt/neutral.top Build_system/. 
        #cp Solvent_npt/*.itp Build_system/. 
        cp Relax_SAM/nvt_relax.nopbc.gro Build_system/. 
        cd Build_system

        zdim=`tail -n1 solvent_npt.gro | awk '{print $3}'`
        ydim=`tail -n1 solvent_npt.gro | awk '{print $2}'`
        xdim=`tail -n1 solvent_npt.gro | awk '{print $1}'`
        zshift=`cat nvt_relax.nopbc.gro | grep LIG | grep C10 | awk '{print $6}' | sort -n | tail -n1`
        zshift=`echo "$zshift * -1" | bc -l`
    
        cp nvt_relax.nopbc.gro bottom_boxed.gro 
        #gmx editconf -f nvt_relax.nopbc.gro \
        #    -translate 0 0 $zshift \
        #    -o bottom_boxed.gro >> $logFile 2>> $errFile 
        #check bottom_boxed.gro 
    
        zdim=`echo "$zdim - $zshift" | bc -l`

        gmx editconf -f solvent_npt.gro \
            -box $xdim $ydim $zdim \
            -o new_box.gro >> $logFile 2>> $errFile 
        check new_box.gro 

        zshift=`cat new_box.gro | grep SOL | grep OW | awk '{print $6}' | sort -n | uniq | tail -n1`
        zshift=`echo "$zdim - $zshift" | bc -l`
    
        gmx editconf -f new_box.gro \
            -translate 0 0 $zshift \
            -o shifted.gro >> $logFile 2>> $errFile 
        check shifted.gro 

        head -n1 shifted.gro > system.gro 
        protAtoms=`cat shifted.gro | tail -n +2 | head -n 1`
        samAtoms=`cat bottom_boxed.gro | tail -n +2 | head -n 1`
        totAtoms=`echo "$protAtoms + $samAtoms" | bc -l`

        echo $totAtoms >> system.gro 
        tail -n +3 shifted.gro | sed '$ d' | grep -v 'SOL' | grep -v 'Cl-' >> system.gro 
        tail -n +3 bottom_boxed.gro | sed '$ d' >> system.gro 
        tail -n +3 shifted.gro | grep 'Cl-' >> system.gro 
        tail -n +3 shifted.gro | grep 'SOL' >> system.gro 
        tail -n1 shifted.gro >> system.gro 

        ## 3, 3 -- None, None for terimini options
        echo '3 3' | gmx pdb2gmx -f system.gro \
            -p neutral.top \
            -ter \
            -ff $forceField \
            -water tip3p \
            -merge interactive \
            -o system.gro >> $logFile 2>> $errFile
        check neutral.top 

        includeLine=`cat -n neutral.top | grep '; Include Position restraint file' | awk '{print $1}'`
        head -n $includeLine neutral.top > system.top  
        echo '#include "distance_restraints.itp"' >> system.top 
        echo '#ifdef POSSULRES' >> system.top 
        echo '#include "posre_SUL.itp"' >> system.top 
        echo '#endif' >> system.top 
        echo >> system.top 
        ((includeLine++)) 
        tail -n +$includeLine neutral.top >> system.top 

        ##Sulfer atoms have changed numbers. Need to rebuild posre_SUL.itp 
        echo "[ position_restraints ]" > posre_SUL.itp 
        echo ";; Pin sulfur atoms to spacing found on SAM" >> posre_SUL.itp 
        for atom in `grep LIG system.gro | grep S1 | awk '{print $3}'` ; do 
            printf "%6i%6i%10.f%10.f%10.f\n" $atom 1 $sulRest $sulRest $sulRest >> posre_SUL.itp 
        done 

        ## We need to rebuild posre.itp. pdb2gmx assigns all heavy atoms to restraints. 
        ## We want the decanethiol heavy atoms to be able to relax with protein position restraints. 
        echo "; Only protein heavy atom position restraints" > posre.itp 
        echo "[ position_restraints ]" >> posre.itp 
        echo "; atom  type      fx      fy      fz" >> posre.itp 
        for atom in `tail -n +2 system.gro | sed '$ d' | grep -v HOH | grep -v CL | grep -v LIG | awk '{print $2"\t"$3}' | grep "^[CNOS]" | awk '{print $2}'` ; do 
            printf "%6i%6i%10.f%10.f%10.f\n" $atom 1 1000 1000 1000 >> posre.itp 
            done 

        add_restraints

        check system.top 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

add_restraints(){ 
    touch empty.ndx 
    echo "r LIG & a C10 
    r GLY & ri 6 & a CA
    r GLY & ri 17 & a CA
    q" > selection.dat 
    
    cat selection.dat | gmx make_ndx -f system.gro \
        -n empty.ndx \
        -o index.ndx >> $logFile 2>> $errFile 
    check index.ndx 
    
    echo '0 1' | gmx mindist -f system.gro \
        -n index.ndx \
        -o gly6_dist.out >> $logFile 2>> errFile 
    check gly6_dist.out 
    
    C10=`cat gly6_dist.out | awk '{print $2}'`
    CA=`cat gly6_dist.out | awk '{print $3}'`

    echo "[ bonds ]" > distance_restraints.itp 
    printf ";%6s%6s%6s%8s%8s%8s%12s\n" ai aj func b0 kb >> distance_restraints.itp
    printf "%6s%6s%6s%8s%8s%8s%12s\n" $C10 $CA 6 $glyDist $glyRest >> distance_restraints.itp

    echo '0 2' | gmx mindist -f system.gro \
        -n index.ndx \
        -o gly17_dist.out >> $logFile 2>> $errFile 
    check gly17_dist.out 

    C10=`cat gly17_dist.out | awk '{print $2}'`
    CA=`cat gly17_dist.out | awk '{print $3}'`
    printf "%6s%6s%6s%8s%8s%8s%12s\n" $C10 $CA 6 $glyDist $glyRest >> distance_restraints.itp
} 

system_steep(){
    printf "\t\tSystem steep.............................." 
    if [ ! -f System_steep/system_steep.gro ] ; then 
        create_dir System_steep
        
        cp Build_system/system.top System_steep/.
        cp Build_system/system.gro System_steep/.
        cp Build_system/*.itp System_steep/. 
        #cp Relax_SAM/*.itp System_steep/. 
        cd System_steep

        gmx grompp -f $MDP/system_steep.mdp \
            -p system.top \
            -c system.gro \
            -o system_steep.tpr >> $logFile 2>> $errFile 
        check system_steep.tpr 

        gmx mdrun -deffnm system_steep >> $logFile 2>> $errFile 
        check system_steep.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

system_nvt(){
    printf "\t\tSystem NVT relaxation....................." 
    if [ ! -f System_nvt/system_nvt.gro ] ; then 
        create_dir System_nvt
        
        cp System_steep/system.top System_nvt/.
        cp System_steep/system_steep.gro System_nvt/.
        cp System_steep/*.itp System_nvt/. 
        cd System_nvt

        gmx grompp -f $MDP/system_nvt_relax.mdp \
            -p system.top \
            -c system_steep.gro \
            -o system_nvt.tpr >> $logFile 2>> $errFile 
        check system_nvt.tpr 

        gmx mdrun -deffnm system_nvt >> $logFile 2>> $errFile 
        check system_nvt.gro 
    
    
        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

production(){
    printf "\t\tProduction run............................" 
    if [ ! -f Production/$MOLEC.nopbc.gro ] ; then 
        create_dir Production
        
        cp System_nvt/system.top Production/.
        cp System_nvt/system_nvt.gro Production/.
        cp System_nvt/*.itp Production/. 
        cd Production

        if [ ! -f $MOLEC.gro ] ; then 
            if [ ! -f $MOLEC.tpr ] ; then 
                gmx grompp -f $MDP/production_sam.mdp \
                    -p system.top \
                    -c system_nvt.gro \
                    -o $MOLEC.tpr >> $logFile 2>> $errFile 
                fi 
                check $MOLEC.tpr 

            if [ -f $MOLEC.cpt ] ; then 
                gmx mdrun -deffnm $MOLEC -cpi $MOLEC.cpt >> $logFile 2>> $errFile  
            else 
                gmx mdrun -deffnm $MOLEC >> $logFile 2>> $errFile 
                fi 
            fi 
        check $MOLEC.gro 

        if [ ! -f $MOLEC.nopbc.xtc ] ; then 
            echo 'System' | gmx trjconv -f $MOLEC.xtc \
                -s $MOLEC.tpr \
                -ur rect \
                -pbc mol \
                -o $MOLEC.nopbc.xtc >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.nopbc.xtc 

        if [ ! -f $MOLEC.nopbc.gro ] ; then 
            echo 'System' | gmx trjconv -f $MOLEC.gro \
                -s $MOLEC.tpr \
                -ur rect \
                -pbc mol \
                -o $MOLEC.nopbc.gro >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.nopbc.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

dssp(){
    printf "\t\tRunning dssp analysis....................." 
    if [ ! -f dssp/helen.nrt ] ; then 
        create_dir dssp
        cd dssp
        clean ##clean early. One of the outputs of gmx do_dssp is a *.dat file. We don't want to delete this while cleaning. 

        #echo 'Protein' | gmx do_dssp -f ../Production/$MOLEC.xtc \
        #    -s ../Production/$MOLEC.tpr \
        #    -ver 1 \
        #    -sss HGI \
        #    -ssdump ssdump.dat \
        #    -o ss.xpm \
        #    -a area.xpm \
        #    -ta totarea.xvg \
        #    -aa averarea.xvg \
        #    -sc scount.xvg >> $logFile 2>> $errFile
        #check scount.xvg ss.xpm area.xpm 

        #gmx xpm2ps -f area.xpm \
        #    -by 10 \
        #    -o area.eps >> $logFile 2>> $errFile 
        #check area.eps

        #gmx xpm2ps -f ss.xpm \
        #    -by 10 \
        #    -o ss.eps >> $logFile 2>> $errFile 
        #check ss.eps

        ##Cut out helen.nrt from scount.xvg. 
        echo "#!/usr/bin/env python
import numpy as np  

with open('scount.xvg') as f : 
    filelines = f.readlines() 

header=0
for line in filelines : 
    if line.startswith('#') or line.startswith('@') : 
        header += 1
        if 'A-Helix' in line : 
            alpha=line.split()[1][1:] 
        if '3-Helix' in line : 
            three=line.split()[1][1:] 
        if '5-Helix' in line : 
            five =line.split()[1][1:] 
header -= 2

data = np.genfromtxt('scount.xvg',skip_header=header)
col1 = data[:,0]

try : 
    col2 = data[:,int(alpha)]
except NameError : 
    col2 = np.zeros(len(data[:,0])) 
try : 
    col3a = data[:,int(three)]
except NameError : 
    col3a = np.zeros(len(data[:,0])) 
try : 
    col3b = data[:,int(five)]
except NameError : 
    col3b = np.zeros(len(data[:,0])) 
col3=col3a + col3b
col4 = 18 - col2 - col3 

for i in range(len(col1)) : 
    print \"%5i%5i%5i%5i\"%(col1[i],col2[i],col3[i],col4[i])" > cut_helen.py 

        python cut_helen.py > helen.nrt 
        check helen.nrt 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

rgyr(){
    printf "\t\tCalculating radius of gyration............" 
    if [ ! -f rgyr/gyrate.xvg ] ; then 
        create_dir rgyr
        cd rgyr
        clean 

        echo 'Protein' | gmx gyrate -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -o gyrate.xvg >> $logFile 2>> $errFile 
        check gyrate.xvg

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

minimage(){
    printf "\t\tCalculating minimum image................." 
    if [ ! -f minimage/mindist.xvg ] ; then 
        create_dir minimage
        cd minimage
        clean 
        
        echo 'Protein' | gmx mindist -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -pi \
            -od mindist.xvg >> $logFile 2>> $errFile 
        check mindist.xvg 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

printf "\n\t\t*** Program Beginning ***\n\n"
cd $MOLEC
build_SAM
layer_relax
protein_steep
solvate
solvent_steep
solvent_nvt
solvent_npt
build_system
system_steep
system_nvt
production
dssp
rgyr 
minimage
cd ../

printf "\n\n\t\t*** Program Ending    ***\n\n"






