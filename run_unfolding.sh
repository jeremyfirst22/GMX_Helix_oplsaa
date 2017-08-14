#!/bin/bash

dim=7.455

usage(){
    echo "USAGE: $0 <PDB file {molec.pdb} > < simulation time (ns) [default=50ns]>"
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
    MOLEC='prep'
    #MOLEC=$(basename $fileName)
    #MOLEC=${MOLEC%.*}_water
else 
    echo "ERROR: Input file must be PDB file (*.pdb)" 
    exit 
    fi 
if [ ! -z $2 ] ; then 
    totSimTime=$2
else 
    totSimTime=50
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

solvate(){
    printf "\t\tSolvating protein........................." 
    if [ ! -f Solvate/neutral.top ] ; then 
        create_dir Solvate
        
        cp Protein_steep/protein_steep.gro Solvate/. 
        cp Protein_steep/$MOLEC.top Solvate/. 
        cd Solvate

        gmx solvate -cp protein_steep.gro \
            -box $dim $dim $dim \
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

heating(){
    printf "\t\tHeating solution to 900K.................." 
    if [ ! -f Heating/heating.gro ] ; then 
        create_dir Heating
        
        cp Solvent_npt/solvent_npt.gro Heating/. 
        cp Solvent_npt/neutral.top Heating/. 
        cp Solvent_npt/*.itp Heating/. 
        cd Heating

        gmx grompp -f $MDP/prep_heating.mdp \
            -c solvent_npt.gro \
            -p neutral.top \
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

heated_nvt(){
    printf "\t\tNVT simulation at 900K...................." 
    if [ ! -f Heated_nvt/heated_nvt.gro ] ; then 
        create_dir Heated_nvt
        
        cp Heating/heating.gro Heated_nvt/. 
        cp Heating/neutral.top Heated_nvt/. 
        cp Heating/*.itp Heated_nvt/. 
        cd Heated_nvt

        gmx grompp -f $MDP/prep_heated.mdp \
            -c heating.gro \
            -p neutral.top \
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
        cp Heated_nvt/neutral.top Cooling/. 
        cp Heated_nvt/*.itp Cooling/. 
        cd Cooling

        gmx grompp -f $MDP/prep_cooling.mdp \
            -c heated_nvt.gro \
            -p neutral.top \
            -o cooling.tpr >> $logFile 2>> $errFile 
        check cooling.tpr

        gmx mdrun -deffnm cooling >> $logFile 2>> $errFile 
        check cooling.gro 

        echo "Protein-H Protein-H" | gmx trjconv -s cooling.tpr \
            -f cooling.gro \
            -o unfolded.pdb \
            -center \
            -ur compact \
            -pbc mol 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

extract(){
    echo "Protein-H Protein-H" | gmx trjconv -s Cooling/cooling.tpr \
        -f Cooling/cooling.gro \
        -o unfolded.pdb \
        -center \
        -ur compact \
        -pbc mol >> $logFile 2>> $errFile 
    check unfolded.pdb 

    printf "\n\t\tUnfolded structure extracted to unfoleded.pdb\n" 
}

printf "\n\t\t*** Program Beginning $MOLEC $totSimTime (ns)***\n\n" 
cd $MOLEC
protein_steep
solvate
solvent_steep
solvent_nvt
solvent_npt
heating
heated_nvt
cooling
extract
cd ../

printf "\n\n\t\t*** Program Ending    ***\n\n" 
