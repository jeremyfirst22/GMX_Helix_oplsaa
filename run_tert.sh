#!/bin/bash
## Jeremy First

TOP=${PWD}
MDP=$TOP/mdp_files
FF=$TOP/GMXFF
forceField=oplsaa
dim=7.600 
fileName=$TOP/StartingStructures/folded.pdb
waterFileName=$TOP/StartingStructures/tip3p.pdb 
tbaFileName=$TOP/StartingStructures/tba.pdb 
molRatio=2
totSimTime=50
SOL=tert
MOLEC=folded_$SOL

verbose=false
analysis=false

usage(){
    echo "USAGE: $0 -f prep/folded/unfolded [ options ]" 
    exit 
}

HELP(){
    echo 
    echo "This program runs molecular dynamic simulations of a peptide in a TBA:Water mixture"
    echo 
    echo "Usage: $0 [options] "
    echo "  -f   folded, unfolded, or prep. Mandatory" 
    echo "  -c   Starting structure  (folded) PDB file: Default = StartingStructures/folded.pdb"
    echo "  -t   Maximum simulation time (ns) : Default = 50 "
    echo "  -a   Perform analyis on trajectory : Default = no"
    echo "  -w   Water structure PDB file: Default = StartingStrcutures/tip3p.pdb" 
    echo "  -b   TBA structure PDB file: Default = StartingStrcutures/tba.pdb" 
    echo "  -R   Volume ratio of TBA:Water in binary solvent: Default = 2"
    echo "  -d   Box dimension (nm) : Default = $dim " 
    echo "  -m   Location of the mdp_files : Default = mdp_files"
    echo "  -p   Location of force field files. : Default = GMXFF"
    echo "  -n   Name of force field : Default = oplsaa"
    echo "  -v   Print all options and quit." 
    echo "  -h   print this usage and exit "
    echo "" 
    exit
}

while getopts :f:c:t:aw:b:R:d:m:p:n:vh opt; do 
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
      w) 
        waterFileName=${TOP}/$OPTARG
        ;; 
      b) 
        tbaFileName=${TOP}/$OPTARG
        ;; 
      R) 
        molRatio=$OPTARG
        ;; 
      d) 
        dim=$OPTARG
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
        fi 
    cd ../
    printf "\n\n\t\t*** Program Ending    ***\n\n" 
}

prep(){
    if [ ! -d prep ] ; then mkdir prep ; fi 
    cd prep
    protein_steep
    prepare_box
    binary_steep
    binary_nvt
    binary_npt
    solvate
    solvent_steep
    solvent_nvt
    solvent_npt
    heating
    heated
    cooling
    cd ../
}

analysis(){
    if [ ! -d $fold ] ; then mkdir $fold ; fi 
    cd $fold 
    dssp
    rgyr
    minimage
#    rdf
    nopbc
    rmsd
    cd_spectra
    cd ../
}

checkInput(){
    if $verbose ; then 
        echo "Folded state : $fold" 
        echo "Input file name: $fileName"
        echo "Water structure name: $waterFileName" 
        echo "TBA structure name: $tbaFileName" 
        echo "Max simultaiton time: $totSimTime"
        echo "Volume ratio of TBA:Water: $molRatio"
        echo "Box dimension : $dim " 
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
    if [ ! -f $waterFileName ] ; then 
        echo "ERROR: $waterFileName not found " 
        exit 
        fi 
    if [[ $waterFileName != *.pdb ]] ; then 
        echo "ERROR: Input file must be PDB file (*.pdb)" 
        exit 
        fi 
    if [ ! -f $tbaFileName ] ; then 
        echo "ERROR: $tbaFileName not found " 
        exit 
        fi 
    if [[ $tbaFileName != *.pdb ]] ; then 
        echo "ERROR: Input file must be PDB file (*.pdb)" 
        exit 
        fi 
    if [[ $molRatio < 0 ]] ; then 
        echo "ERROR: Box dimension must be larger than 0"   
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
    if [[ $dim < 0 ]] ; then 
        echo "ERROR: Box dimension must be larger than 0" 
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
        
        cp $waterFileName Prepare_box/tip3p.pdb 
        cp $tbaFileName Prepare_box/tba.pdb 
        cp Protein_steep/$MOLEC.top Prepare_box/. 
        cp Protein_steep/protein_steep.gro Prepare_box/. 
        cd Prepare_box/. 
        check tip3p.pdb tba.pdb 

        ### 0.800 g/cm^3, density of 2:1 TBA:H20
        ### 166.255 g/mol, molar mass of one mol of 2 tba molecules and 1 water molecule
        #numWat=`echo "print int(round(($dim * 10 **-7) ** 3 * 0.800 / 166.255 * 6.022 * 10 ** 23))" | python`
        #numTBA=`echo "$numWat * 2 " | bc -l`
        
        echo "#!/usr/bin/env python
Na = 6.022 * 10 ** 23 #Avogadro's number 
pTBA = 0.77014 ## g/cm**3 
mmTBA = 74.12 ## g/mol 
pH2O = 0.994029 ## g/cm**3 
mmH2O = 18.01528 ## g/mol
molRatio = $molRatio ##vol Ratio of tba:water
V = ($dim * 10 ** -7)**3 ##cm**3
xTBA = molRatio / (molRatio + 1.) ##vol fraction
xH2O = 1 - xTBA ##vol fraction
nTBA = xTBA * V * pTBA * Na / mmTBA 
nH2O = xH2O * V * pH2O * Na / mmH2O 
print \"%i\t%i\"%(round(nTBA),round(nH2O))" > calc_num_molecules.py 
        python calc_num_molecules.py > molecules.xvg 
        numTBA=`awk '{print $1}' molecules.xvg`
        numWat=`awk '{print $2}' molecules.xvg`

        gmx insert-molecules -ci tba.pdb \
            -box $dim $dim $dim \
            -allpair \
            -nmol $numTBA \
            -o tert_box.gro >> $logFile 2>> $errFile 
        check tert_box.gro 

        gmx insert-molecules -ci tip3p.pdb \
            -f tert_box.gro \
            -allpair \
            -nmol $numWat \
            -o mixture.gro >> $logFile 2>> $errFile 
        check mixture.gro 

        gmx pdb2gmx -f mixture.gro \
            -ff $forceField \
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

        gmx editconf -f solvated.gro \
            -o solvated.pdb >> $logFile 2>> $errFile 
        check solvated.pdb 

        awk '/H2  NH2/{print;print "TER";next}1' solvated.pdb > temp.pdb
        mv temp.pdb solvated.pdb

        ## 3, 3 -- None, None for terimini options
        echo '3 3' | gmx pdb2gmx -f solvated.pdb \
            -ff $forceField \
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
        
        echo 'Water' | gmx genion -s genion.tpr \
            -neutral \
            -nname 'CL' \
            -pname 'NA' \
            -o neutral.gro >> $logFile 2>> $errFile 
        check neutral.gro 

        gmx editconf -f neutral.gro \
            -o neutral.pdb >> $logFile 2>> $errFile
        check neutral.pdb

        awk '/H2  NH2/{print;print "TER";next}1' neutral.pdb > temp.pdb
        mv temp.pdb neutral.pdb 

        ## 3, 3 -- None, None for terimini options
        echo '3 3' | gmx pdb2gmx -f neutral.pdb \
            -ff $forceField \
            -water tip3p \
            -p neutral.top \
            -ter \
            -o neutral.gro >> $logFile 2>> $errFile 
        check neutral.top 

        ##Need to update so that TBUT and ions don't get restrained with protein
        sed 's/POSRES/POSRES_SOLVENT/' neutral_Solvent2.itp > temp.itp 
        mv temp.itp neutral_Solvent2.itp 
        check neutral_Solvent2.itp 

        sed 's/POSRES/POSRES_IONS/' neutral_Ion3.itp > temp.itp 
        mv temp.itp neutral_Ion3.itp
        check neutral_Ion3.itp 

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

heated(){
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
    if [[ ! -s prep/Cooling/cooling.gro || ! -s prep/Solvent_npt/solvent_npt.gro ]] ; then 
        echo ; echo 
        echo "ERROR: You must run prep first!" 
        fi 
    check prep/Cooling/cooling.gro prep/Solvent_npt/solvent_npt.gro 

    if [ ! -d $fold ] ; then mkdir $fold ; fi 
    cd $fold 
    MOLEC=${fold}_${SOL}

    if [ ! -f Production/${MOLEC}_${totSimTime}ns.gro ] ; then 
        printf "\n" 
        create_dir Production
        
        cp ../prep/Solvent_npt/neutral.top Production/.
        cp ../prep/Solvent_npt/*.itp Production/. 

        if [ $fold = "folded" ] ; then 
            cp ../prep/Solvent_npt/solvent_npt.gro Production/.
            startStructure=solvent_npt.gro
            MOLEC=folded_$SOL
        else 
            cp ../prep/Cooling/cooling.gro Production/.
            startStructure=cooling.gro 
            MOLEC=unfolded_$SOL
            fi 

        cd Production

        if [ ! -f $MOLEC.tpr ] ; then 
            gmx grompp -f $MDP/production_solvent.mdp \
                -p neutral.top \
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
    cd ../
} 

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
    col2 = data[:,int(alpha)+1]
except NameError : 
    col2 = np.zeros(len(data[:,0])) 
try : 
    col3a = data[:,int(three)+1]
except NameError : 
    col3a = np.zeros(len(data[:,0])) 
try : 
    col3b = data[:,int(five)+1]
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
        echo "r TBUT & a C" > selection.dat 
        echo "r SOL & a OW" >> selection.dat 
        echo "a CD1 or a CD2 & r LEU" >> selection.dat 
        echo "r LYS & a NZ" >> selection.dat 
        echo "q" >> selection.dat 

        cat selection.dat | gmx make_ndx -f ../Production/solvent_npt.gro \
            -n empty.ndx \
            -o index.ndx >> $logFile 2>> $errFile 
        check index.ndx 

        echo '0 0' | gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -o tba_tba.xvg >> $logFile 2>> $errFile 
        check tba_tba.xvg

        echo '0 1' | gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -o tba_wat.xvg >> $logFile 2>> $errFile 
        check tba_wat.xvg 

        echo '1 1' | gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -o wat_wat.xvg >> $logFile 2>> $errFile 
        check wat_wat.xvg

        echo '2 0' | gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -o leu_tba.xvg >> $logFile 2>> $errFile 
        check leu_tba.xvg

        echo '2 1' | gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -o leu_wat.xvg >> $logFile 2>> $errFile 
        check leu_wat.xvg 

        echo '3 0' | gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -o lys_tba.xvg >> $logFile 2>> $errFile 
        check lys_tba.xvg 

        echo '3 1' | gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -o lys_wat.xvg >> $logFile 2>> $errFile 
        check lys_wat.xvg 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

nopbc(){
    printf "\t\tCreating smoothed trajectory.............." 
    if [ ! -f nopbc/smooth.xtc ] ; then 
        create_dir nopbc
        cd nopbc
        clean 

        if [ $fold = "folded" ] ; then 
            startStructure=solvent_npt.gro 
            MOLEC=folded_$SOL
        else 
            startStructure=cooling.gro 
            MOLEC=unfolded_$SOL
            fi 

        echo "Protein Protein" | gmx trjconv -f ../Production/$startStructure \
            -s ../Production/$MOLEC.tpr \
            -pbc mol \
            -ur compact \
            -center \
            -o nopbc.gro >> $logFile 2>> $errFile 
        check nopbc.gro 

        echo "Protein Protein" | gmx trjconv -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -pbc mol \
            -ur compact \
            -center \
            -dt 40 \
            -o nopbc.xtc >> $logFile 2>> $errFile 
        check nopbc.xtc 

        echo "Backbone Protein" | gmx trjconv -f nopbc.xtc \
            -s ../Production/$MOLEC.tpr \
            -fit rot+trans \
            -o fit.xtc >> $logFile 2>> $errFile 
        check fit.xtc

        echo "Backbone Protein" | gmx trjconv -f fit.xtc \
            -s ../Production/$MOLEC.tpr \
            -fit progressive \
            -o smooth.xtc >> $logFile 2>> $errFile 
        check smooth.xtc 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

rmsd(){
    printf "\t\tCalculating RMSD.........................." 
    if [ ! -f rmsd/rmsd.xvg ] ; then 
        create_dir rmsd
        cd rmsd
        clean 

        echo "Backbone Backbone" | gmx rms -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -o rmsd.xvg >> $logFile 2>> $errFile 
        check rmsd.xvg  

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

cd_spectra(){
    printf "\t\tExtracting frames for DichroCalc.........." 
    if [ ! -f cd_spectra/extract.tar.gz ] ; then 
        create_dir cd_spectra
        cd cd_spectra
        clean 
        create_dir extract
        cd extract
        clean 

        echo 'Protein' | gmx trjconv -f ../../Production/$MOLEC.xtc \
            -s ../../Production/$MOLEC.tpr \
            -dt 100 \
            -pbc mol \
            -sep \
            -o frame.pdb >> $logFile 2>> $errFile 
        cd ../

        tar cvfz extract.tar.gz extract >> $logFile 2>> $errFile 
        check extract.tar.gz 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

main
