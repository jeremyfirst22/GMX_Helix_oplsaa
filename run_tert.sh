#!/bin/bash

dim=5.964

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
    MOLEC=${MOLEC%.*}_tert
else 
    echo "ERROR: Input file must be PDB file (*.pdb)" 
    exit 
    fi 
if [ ! -d $MOLEC ] ; then mkdir $MOLEC ; fi 
if [ ! -f $MOLEC/$fileName ] ; then cp $fileName $MOLEC/$MOLEC.pdb ; fi 

if [[ ! -f StartingStructures/tip3p.pdb || ! -f StartingStructures/tba.pdb ]] ; then 
    echo "ERROR: tip3p.pdb and tba.pdb required in StartingStrcutures for this solvent" 
    exit 
    fi 
cp StartingStructures/tip3p.pdb $MOLEC/.
cp StartingStructures/tba.pdb $MOLEC/.

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

check(){
   for arg in $@ ; do 
        if [ ! -s $arg ] ; then 
            echo ; echo "ERROR: $arg missing. Exitting" 
            exit 
            fi 
        done 
}

clean(){
    if [ -d oplsaa.ff ] ; then rm -r oplsaa.ff *.dat ; fi 
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

        xdim=$dim
        ydim=$dim
        zdim=$dim
            #-d 1.0 \   ## Box size is bigger than -d 1.0 nm, and consistent with sam_box.
        echo 'Backbone' | gmx editconf -f $MOLEC.gro \
            -box $xdim $ydim $zdim \
            -bt tric \
            -o boxed.gro >> $logFile 2>> $errFile
        check boxed.gro 

        gmx grompp -f $MDP/protein_steep.mdp \
            -c boxed.gro \
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

prepare_box(){
    printf "\t\tPreparing binary mixture.................." 
    if [ ! -f Prepare_box/mixture.top ] ; then 
        create_dir Prepare_box
        
        cp tip3p.pdb Prepare_box/. 
        cp tba.pdb Prepare_box/. 
        cp Protein_steep/$MOLEC.top Prepare_box/. 
        cp Protein_steep/protein_steep.gro Prepare_box/. 
        cd Prepare_box/. 
        check tip3p.pdb tba.pdb 

        gmx insert-molecules -ci tba.pdb \
            -box $dim $dim $dim \
            -nmol 1230 \
            -o tert_box.gro >> $logFile 2>> $errFile 
        check tert_box.gro 

        gmx insert-molecules -ci tip3p.pdb \
            -f tert_box.gro \
            -nmol 615 \
            -o mixture.gro >> $logFile 2>> $errFile 
        check mixture.gro 

        gmx pdb2gmx -f mixture.gro \
            -ff oplsaa \
            -water tip3p \
            -o mixture.gro \
            -p mixture.top >> $logFile 2>> $errFile 
        check mixture.top 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

binary_steep(){
    printf "\t\tBinary mixture steep......................" 
    if [ ! -f Binary_steep/binary_steep.gro ] ; then 
        create_dir Binary_steep
        
        cp Prepare_box/mixture.gro Binary_steep/. 
        cp Prepare_box/mixture.top Binary_steep/. 
        cd Binary_steep

        gmx grompp -f $MDP/binary_steep.mdp \
            -p mixture.top \
            -c mixture.gro \
            -o binary_steep.tpr >> $logFile 2>> $errFile  
        check binary_steep.tpr 

        gmx mdrun -deffnm binary_steep >> $logFile 2>> $errFile 
        check binary_steep.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

binary_nvt(){
    printf "\t\tBinary mixture NVT relaxation............." 
    if [ ! -f Binary_nvt/binary_nvt.gro ] ; then 
        create_dir Binary_nvt
        
        cp Binary_steep/binary_steep.gro Binary_nvt/. 
        cp Binary_steep/mixture.top Binary_nvt/. 
        cd Binary_nvt

        gmx grompp -f $MDP/binary_nvt_relax.mdp \
            -p mixture.top \
            -c binary_steep.gro \
            -o binary_nvt.tpr >> $logFile 2>> $errFile  
        check binary_nvt.tpr 

        gmx mdrun -deffnm binary_nvt >> $logFile 2>> $errFile 
        check binary_nvt.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

binary_npt(){
    printf "\t\tBinary mixture NPT relaxation............." 
    if [ ! -f Binary_npt/binary_npt.gro ] ; then 
        create_dir Binary_npt
        
        cp Binary_nvt/binary_nvt.gro Binary_npt/. 
        cp Binary_nvt/mixture.top Binary_npt/. 
        cd Binary_npt

        gmx grompp -f $MDP/binary_npt_relax.mdp \
            -p mixture.top \
            -c binary_nvt.gro \
            -o binary_npt.tpr >> $logFile 2>> $errFile 
        check binary_npt.tpr

        gmx mdrun -deffnm binary_npt >> $logFile 2>> $errFile 
        check binary_npt.gro 

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
        cp Binary_npt/binary_npt.gro Solvate/. 
        cd Solvate

        gmx solvate -cp protein_steep.gro \
            -cs binary_npt.gro \
            -o solvated.gro >> $logFile 2>> $errFile 
        check solvated.gro

        ## 3, 3 -- None, None for terimini options
        echo '3 3' | gmx pdb2gmx -f solvated.gro \
            -ff oplsaa \
            -water tip3p \
            -p solvated.top \
            -ter \
            -o solvated.gro >> $logFile 2>> $errFile 
        check solvated.top 

        gmx grompp -f $MDP/vac_md.mdp \
            -p solvated.top \
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
            -ff oplsaa \
            -water tip3p \
            -p neutral.top \
            -ter \
            -o neutral.gro >> $logFile 2>> $errFile 
        check neutral.top 

        ##Need to rebuild posre_Protein.itp to be only heavy atoms from protein
        ## pdb2gmx assigns TBUT heavy atoms to restraints called by POSRES, but 
        ##    we want all solvent molecuels to relax. 
        rm posre.itp 
        echo "; Position restraints ONLY for protein heavy atoms" > posre_Protein.itp 
        echo "[ position_restraints ]" >> posre_Protein.itp 
        echo "; atom  type      fx      fy      fz" >> posre_Protein.itp 
        for atom in `tail -n +3 neutral.gro | sed '$ d' | grep -v HOH | grep -v TBUT | grep -v CL | awk '{print $2"\t"$3}' | grep "^[CNOS]" | awk '{print $2}'` ; do 
            printf "%6i%6i%10.f%10.f%10.f\n" $atom 1 1000 1000 1000 >> posre_Protein.itp 
            done 

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
    printf "\t\tSolvent NPT isotropic relaxation.........." 
    if [ ! -f Solvent_npt/solvent_npt.gro ] ; then 
        create_dir Solvent_npt
        
        cp Solvent_nvt/solvent_nvt.gro Solvent_npt/. 
        cp Solvent_nvt/neutral.top Solvent_npt/. 
        cp Solvent_nvt/*.itp Solvent_npt/. 
        cd Solvent_npt

        gmx grompp -f $MDP/solvent_isotropic_npt_relax.mdp \
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

production(){
    printf "\t\tProduction run............................" 
    if [ ! -f Production/$MOLEC.nopbc.gro ] ; then 
#            CA1=$(grep " CA " $MOLEC.npt_relax.gro | grep "GLY" | awk '{print $3}' | head -n1) 
#            CA2=$(grep " CA " $MOLEC.npt_relax.gro | grep "GLY" | awk '{print $3}' | tail -n1) 
#            if [[ -z $CA1 || -z $CA2 ]] ; then echo "ERROR: Failed to find CA1 or CA2" ; exit ; fi 

#            echo "[ bonds ]" > distance_restraints.itp 
#            printf ";%6s%6s%6s%8s%8s%8s%12s\n" ai aj func b0 kb >> distance_restraints.itp 
#            printf "%6s%6s%6s%8s%8s%8s%12s\n" $CA1 $CA2 6 1.8 1000 >> distance_restraints.itp 

#            if ! grep -sq "distance_restraints.itp" $MOLEC.restraint.top ; then 
#                awk -v molec=$MOLEC '/.neutral_Protein.itp\"/{print;print"#include \"distance_restraints.itp\" ";next}1' $MOLEC.neutral.top > $MOLEC.restraint.top 
#                fi 

#            if ! grep -sq "distance_restraints.itp" $MOLEC.restraint.top ; then 
#                echo "ERROR: distance restraint not added to topology." 
#                exit 
#                fi 
        create_dir Production
        
        cp Solvent_npt/neutral.top Production/.
        cp Solvent_npt/solvent_npt.gro Production/.
        cp Solvent_npt/*.itp Production/. 
        cd Production

        if [ ! -f $MOLEC.gro ] ; then 
            if [ ! -f $MOLEC.tpr ] ; then 
                gmx grompp -f $MDP/production_sam.mdp \
                    -p neutral.top \
                    -c solvent_npt.gro \
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
            echo 'Protein System' | gmx trjconv -f $MOLEC.xtc \
                -center \
                -s $MOLEC.tpr \
                -ur rect \
                -pbc mol \
                -o $MOLEC.nopbc.xtc >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.nopbc.xtc 

        if [ ! -f $MOLEC.nopbc.gro ] ; then 
            echo 'Protein System' | gmx trjconv -f $MOLEC.gro \
                -center \
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

printf "\n\t\t*** Program Beginning ***\n\n" 
cd $MOLEC
protein_steep
prepare_box
binary_steep
binary_nvt
binary_npt
solvate
solvent_steep
solvent_nvt
solvent_npt
production 
cd ../

printf "\n\n\t\t*** Program Ending    ***\n\n" 
