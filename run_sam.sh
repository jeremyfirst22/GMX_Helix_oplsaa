#!/bin/bash
## Jeremy First

TOP=${PWD}
MDP=$TOP/mdp_files
FF=$TOP/GMXFF
forceField=oplsaa
numMols=14
decanethiolFile=$TOP/StartingStructures/decanethiol.pdb 
spacing=0.497 ##nm 
glyDist=0.6822 ##nm from geometry optimization, Spartan '16, DFT wB97X-D, 6-311+g**
glyRest=1000   #kJ/mol/nm
sulRest=200000 #kJ/mol/nm
fileName=$TOP/StartingStructures/folded.pdb
totSimTime=50
SOL=sam
MOLEC=folded_$SOL

verbose=false
analysis=false

usage(){
    echo "USAGE: $0 -f prep/folded/unfolded [ options ]" 
    exit 
}

HELP(){
    echo 
    echo "This program runs molecular dynamic simulations of a peptide on a SAM surface" 
    echo 
    echo "Usage: $0 [options] "
    echo "  -f   folded, unfolded, or prep. Mandatory" 
    echo "  -c   Starting structure  (folded) PDB file: Default = StartingStructures/folded.pdb"
    echo "  -t   Maximum simulation time (ns) : Default = 50 "
    echo "  -a   Perform analyis on trajectory : Default = no"
    echo "  -n   Number of decanethiol molecules in each direction of the layer : Default = 14" 
    echo "  -d   Decanethiol structure PDB file : Default = StartingStructures/decanethiol.pdb" 
    echo "  -D   Distance (nm) between each decanethiol molecule : Default = 0.497" 
    echo "  -g   Distnace (nm) between the alpha carbon of glycines to C10 of the SAM layer" 
    echo "              Default = 0.6822 from Geoemetry optimization, Spartan '16, DFT wB97X-D, 6-311+g**"
    echo "  -r   Restraint force constant (kJ/mol/nm) of the "bond" from glycine to nearest decanethiol C10"
    echo "              Default = 1000 "
    echo "  -R   Restraint force constant (kJ/mol/nm) of the "bond" keeping sulfur atoms in plane" 
    echo "              Default = 200000" 
    echo "  -m   Location of the mdp_files : Default = mdp_files"
    echo "  -p   Location of force field files. : Default = GMXFF"
    echo "  -n   Name of force field : Default = oplsaa"
    echo "  -v   Print all options and quit." 
    echo "  -h   print this usage and exit "
    echo "" 
    exit
}

while getopts :f:c:t:an:d:D:g:r:R:m:p:n:vh opt; do 
   case $opt in 
      f) 
        fold=$OPTARG
        ;; 
      c)
        fileName=${TOP}/$OPTARG
        ;; 
      t)
        totSimTime=$OPTARG
        ;; 
      a) 
        analysis=true
        ;; 
      d) 
        decanethiolFile=${TOP}/$OPTARG
        ;; 
      D) 
        spacing=$OPTARG
        ;; 
      g) 
        glyDist=$OPTARG
        ;; 
      r) 
        glyRest=$OPTARG
        ;; 
      R) 
        sulRest=$OPTARG
        ;; 
      m) 
        MDP=${TOP}/$OPTARG
        ;; 
      p) 
        FF=${TOP}/$OPTARG
        ;; 
      n) 
        forceField=$OPTARG
        ;; 
      v) 
        verbose=true
        ;; 
      :) 
        echo " option -$OPTARG requires an argument! "
        usage
        ;; 
      h)
        HELP
        ;; 
      \?) 
        echo "Invalid option -$OPTARG "
        HELP 
        ;; 
   esac
   done 

main(){
    logFile=$TOP/$SOL/$fold/$fold.log
    errFile=$TOP/$SOL/$fold/$fold.err
    checkInput
    printf "\n\t\t*** Program Beginning $SOL_$fold $totSimTime (ns)***\n\n" 
    if [ ! -d $SOL ] ; then mkdir $SOL ; fi 
    cd $SOL
    if [ $fold == "prep" ] ; then 
        prep
    else  
        production 
        if $analysis ; then 
            analysis
            fi 
        cd ../
        fi 
    cd ../
    printf "\n\n\t\t*** Program Ending    ***\n\n" 
}

prep(){
    if [ ! -d prep ] ; then mkdir prep ; fi 
    cd prep
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
    heating
    heated
    cooling
    cd ../
}

analysis(){
    if [ ! -d $fold ] ; then mkdir $fold ; fi 
    dssp
    rgyr
    minimage
    rdf
}

checkInput(){
    if $verbose ; then 
        echo "Folded state : $fold" 
        echo "Input file name: $fileName"
        echo "Max simultaiton time: $totSimTime"
        echo "Number of decanethiols: $numMols" 
        echo "Decanethiol starting structure : $decanethiolFile" 
        echo "Distance between each decanethiol : $spacing"
        echo "Distance of glycine to SAM : $glyDist" 
        echo "Glycine restraint : $glyRest"
        echo "Sulfur restraint : $sulRest" 
        echo "Perform analysis : $analysis " 
        echo "Verbose = $verbose"
        echo "mdp files = $MDP " 
        echo "force field directory = $FF" 
        echo "force field name = $forceField " 
        echo "Log file = $logFile " 
        echo "Error file = $errFile " 
        echo ""
        exit 
        fi 
    
    if [ ! -f $fileName ] ; then 
        echo "ERROR: $fileName not found " 
        exit 
        fi 
    if [[ $fileName != *.pdb ]] ; then 
        echo "ERROR: Input file must be PDB file (*.pdb)" 
        exit 
        fi 
    if [ ! -f $decanethiolFile ] ; then 
        echo "ERROR: $decanethiolFile not found " 
        exit 
        fi 
    if [[ $decanethiolFile != *.pdb ]] ; then 
        echo "ERROR: Input file must be PDB file (*.pdb)" 
        exit 
        fi 
    if [ ! -z $fold ] ; then 
        if [[ $fold != 'folded' && $fold != 'unfolded' && $fold != 'prep' ]] ; then 
            echo "ERROR: $fold not recognized"
            HELP
            fi 
    else 
        echo "ERROR: Must specify folded/unfoled" 
        usage
        fi 
    if [ $((numMols%2)) -ne 0 ] ; then 
        echo "ERROR: numMols must be even for sam layer to tesselate" 
        exit 
        fi 
    if [[ $glyDist < 0 ]] ; then 
        echo "ERROR: Glycine distance to SAM must be larger than zero" 
        exit 
        fi 
    if [[ $spacing < 0 ]] ; then 
        echo "ERROR: spacing between dcanethiols must be larger than zero" 
        exit 
        fi 
    if [[ $glyRest < 0 ]] ; then 
        echo "ERROR: glycine restraint force constant must be larger than zero" 
        exit 
        fi 
    if [[ $sulRest < 0 ]] ; then 
        echo "ERROR: sulfur restraint force constant must be larger than zero" 
        exit 
        fi 
    if [ $totSimTime -le 0 ] ; then 
        echo "ERROR: Total simulation time must be larger than 0" 
        exit 
        fi 
    if [ ! -d $FF/$forceField.ff ] ; then 
        echo ; echo "ERROR: FF not found" 
        exit
        fi  
    if [ ! -d $MDP ] ; then 
        echo ; echo "ERROR: mdp files directory not found" 
        exit 
        fi 
    check $MDP $fileName $FF $FF/$forceField.ff 
} 

check(){
   for arg in $@ ; do 
        if [ ! -s $arg ] ; then 
            echo ; echo "ERROR: $arg missing. Exitting" 
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
            cp $FF/*.dat $dirName/. 
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
        
        #cp $MOLEC.pdb Build_SAM/. 
        cp $decanethiolFile Build_SAM/decanethiol.pdb 
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
            printf "%6i%6i%10.f%10.f%10.f\n" $atom 1 0 0 $sulRest >> posre_SUL.itp 
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

        echo 'System' | gmx trjconv -f nvt_relax.gro \
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
        
        cp $fileName Protein_steep/.
        cd Protein_steep

        ## 3, 3 -- None, None for termini options
        echo '3 3' | gmx pdb2gmx -f $(basename $fileName) \
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
        gmx editconf -f aligned.gro \
            -rotate -60 0 90 \
            -o rotated.gro >> $logFile 2>> $errFile 
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
        
        echo 'Water' | gmx genion -s genion.tpr \
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

        sed 's/POSRES/POSRES_IONS/' neutral_Ion2.itp > temp.itp 
        mv temp.itp neutral_Ion2.itp 

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
        cp Solvate/neutral_*.itp Solvent_steep/. 
        cp Solvate/posre_*.itp Solvent_steep/. 
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

heating(){
    printf "\t\tHeating solution to 900K.................." 
    if [ ! -f Heating/heating.gro ] ; then 
        create_dir Heating
        
        cp System_nvt/system_nvt.gro Heating/. 
        cp System_nvt/system.top Heating/. 
        cp System_nvt/*.itp Heating/. 
        cd Heating

        gmx grompp -f $MDP/prep_heating.mdp \
            -c system_nvt.gro \
            -p system.top \
            -o heating.tpr >> $logFile 2>> $errFile  
        check heating.tpr 

        gmx mdrun -deffnm heating >> $logFile 2>> $errFile 
        check heating.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

heated(){
    printf "\t\tNVT simulation at 900K...................." 
    if [ ! -f Heated_nvt/heated_nvt.gro ] ; then 
        create_dir Heated_nvt
        
        cp Heating/heating.gro Heated_nvt/. 
        cp Heating/system.top Heated_nvt/. 
        cp Heating/*.itp Heated_nvt/. 
        cd Heated_nvt

        gmx grompp -f $MDP/prep_heated.mdp \
            -c heating.gro \
            -p system.top \
            -o heated_nvt.tpr >> $logFile 2>> $errFile 
        check heated_nvt.tpr

        gmx mdrun -deffnm heated_nvt >> $logFile 2>> $errFile 
        check heated_nvt.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

cooling(){
    printf "\t\tCooling solution to 300K.................." 
    if [ ! -f Cooling/cooling.gro ] ; then 
        create_dir Cooling
        
        cp Heated_nvt/heated_nvt.gro Cooling/. 
        cp Heated_nvt/system.top Cooling/. 
        cp Heated_nvt/*.itp Cooling/. 
        cd Cooling

        gmx grompp -f $MDP/prep_cooling.mdp \
            -c heated_nvt.gro \
            -p system.top \
            -o cooling.tpr >> $logFile 2>> $errFile 
        check cooling.tpr

        gmx mdrun -deffnm cooling >> $logFile 2>> $errFile 
        check cooling.gro 

        echo "Protein-H Protein-H" | gmx trjconv -s cooling.tpr \
            -f cooling.gro \
            -o unfolded.pdb \
            -center \
            -ur compact \
            -pbc mol >> $logFile 2>> $errFile

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

production(){
    printf "\t\tProduction run............................" 
    if [[ ! -s prep/Cooling/cooling.gro || ! -s prep/System_nvt/system_nvt.gro ]] ; then 
        echo ; echo 
        echo "ERROR: You must run prep first!" 
        fi 
    check prep/Cooling/cooling.gro prep/System_nvt/system_nvt.gro 

    if [ ! -d $fold ] ; then mkdir $fold ; fi 
    cd $fold 
    MOLEC=${fold}_${SOL}

    if [ ! -f Production/${MOLEC}_${totSimTime}ns.gro ] ; then 
        printf "\n" 
        create_dir Production
        
        cp ../prep/System_nvt/system.top Production/.
        cp ../prep/System_nvt/*.itp Production/. 

        if [ $fold = "folded" ] ; then 
            cp ../prep/System_nvt/system_nvt.gro Production/.
            startStructure=system_nvt.gro 
            MOLEC=folded_$SOL
        else 
            cp ../prep/Cooling/cooling.gro Production/.
            startStructure=cooling.gro 
            MOLEC=unfolded_$SOL
            fi 

        cd Production

        if [ ! -f $MOLEC.tpr ] ; then 
            gmx grompp -f $MDP/production_sam.mdp \
                -p system.top \
            -c $startStructure \
            -o $MOLEC.tpr >> $logFile 2>> $errFile 
        fi 
        check $MOLEC.tpr 

        simTime=0
        while [ $simTime -lt $totSimTime ] ; do 
            ((simTime+=50))
            printf "\t\t\t%10i ns....................." $simTime

            if [ ! -f ${MOLEC}_${simTime}ns.gro ] ; then 
                if [ ! -f $simTime.tpr ] ; then 
                    gmx convert-tpr -s $MOLEC.tpr \
                        -until $((simTime*1000)) \
                        -o $simTime.tpr >> $logFile 2>> $errFile 
                    fi 
                check $simTime.tpr 

                if [ -f $MOLEC.cpt ] ; then 
                    gmx mdrun -deffnm $MOLEC \
                        -s $simTime.tpr \
                        -cpi $MOLEC.cpt >> $logFile 2>> $errFile  
                else 
                    gmx mdrun -deffnm $MOLEC \
                        -s $simTime.tpr >> $logFile 2>> $errFile
                    fi 
                check $MOLEC.gro 
                mv $MOLEC.gro ${MOLEC}_${simTime}ns.gro 
                check ${MOLEC}_${simTime}ns.gro 
                printf "Success\n" 
            else       
                printf "Skipped\n" 
                fi 
            done 

        printf "\n" 
        clean
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

#production(){
#    printf "\t\tProduction run............................" 
#    if [ ! -f Production/${MOLEC}_${totSimTime}ns.gro ] ; then 
#        printf "\n" 
#        create_dir Production
#        
#        cp Solvent_npt/neutral.top Production/.
#        cp Solvent_npt/solvent_npt.gro Production/.
#        cp Solvent_npt/*.itp Production/. 
#        cd Production
#
#        if [ ! -f $MOLEC.tpr ] ; then 
#            gmx grompp -f $MDP/production_solvent.mdp \
#                -p neutral.top \
#            -c solvent_npt.gro \
#            -o $MOLEC.tpr >> $logFile 2>> $errFile 
#        fi 
#        check $MOLEC.tpr 
#
#        simTime=0
#        while [ $simTime -lt $totSimTime ] ; do 
#            ((simTime+=50))
#            printf "\t\t\t%10i ns....................." $simTime
#
#            if [ ! -f ${MOLEC}_${simTime}ns.gro ] ; then 
#                if [ ! -f $simTime.tpr ] ; then 
#                    gmx convert-tpr -s $MOLEC.tpr \
#                        -until $((simTime*1000)) \
#                        -o $simTime.tpr >> $logFile 2>> $errFile 
#                    fi 
#                check $simTime.tpr 
#
#                if [ -f $MOLEC.cpt ] ; then 
#                    gmx mdrun -deffnm $MOLEC \
#                        -s $simTime.tpr \
#                        -cpi $MOLEC.cpt >> $logFile 2>> $errFile  
#                else 
#                    gmx mdrun -deffnm $MOLEC \
#                        -s $simTime.tpr >> $logFile 2>> $errFile
#                    fi 
#                check $MOLEC.gro 
#                mv $MOLEC.gro ${MOLEC}_${simTime}ns.gro 
#                check ${MOLEC}_${simTime}ns.gro 
#                printf "Success\n" 
#            else       
#                printf "Skipped\n" 
#                fi 
#            done 
#
#        printf "\n" 
#        clean
#        cd ../
#    else
#        printf "Skipped\n"
#        fi  
#} 

dssp(){
    printf "\t\tRunning dssp analysis....................." 
    if [ ! -f dssp/helen.nrt ] ; then 
        create_dir dssp
        cd dssp
        clean ##clean early. One of the outputs of gmx do_dssp is a *.dat file. We don't want to delete this while cleaning. 

        echo 'Protein' | gmx do_dssp -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -ver 1 \
            -sss HGI \
            -ssdump ssdump.dat \
            -o ss.xpm \
            -a area.xpm \
            -ta totarea.xvg \
            -aa averarea.xvg \
            -sc scount.xvg >> $logFile 2>> $errFile
        check scount.xvg ss.xpm area.xpm 

        gmx xpm2ps -f area.xpm \
            -by 10 \
            -o area.eps >> $logFile 2>> $errFile 
        check area.eps

        gmx xpm2ps -f ss.xpm \
            -by 10 \
            -o ss.eps >> $logFile 2>> $errFile 
        check ss.eps

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

rdf(){
    printf "\t\tCalculating RDFs.........................." 
    if [ ! -f rdf/lys_wat.xvg ] ; then 
        create_dir rdf
        cd rdf
        clean 

        touch empty.ndx 
        echo "r SOL & a OW" > selection.dat 
        echo "a CD1 or CD2 & r LEU" >> selection.dat 
        echo "r LYS & a NZ" >> selection.dat 
        echo "r LIG & a C10" >> selection.dat 
        echo "q" >> selection.dat 

        cat selection.dat | gmx make_ndx -f ../Production/system_nvt.gro \
            -n empty.ndx \
            -o index.ndx >> $logFile 2>> $errFile 
        check index.ndx 

        echo '0 0' | gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -o wat_wat.xvg >> $logFile 2>> $errFile 
        check wat_wat.xvg

        echo '1 0' | gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -o leu_wat.xvg >> $logFile 2>> $errFile 
        check leu_wat.xvg 

        echo '2 0' | gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -o lys_wat.xvg >> $logFile 2>> $errFile 
        check lys_wat.xvg 

        echo '3 1' | gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -o sam_leu.xvg >> $logFile 2>> $errFile 
        check sam_leu.xvg 

        echo '3 2' | gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -o sam_lys.xvg >> $logFile 2>> $errFile 
        check sam_lys.xvg 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

main
