#!/bin/bash

dim=7.455

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

        ### 0.800 g/cm^3, density of 2:1 TBA:H20
        ### 166.255 g/mol, molar mass of one mol of 2 tba molecules and 1 water molecule
        #numWat=`echo "print int(round(($dim * 10 **-7) ** 3 * 0.800 / 166.255 * 6.022 * 10 ** 23))" | python`
        #numTBA=`echo "$numWat * 2 " | bc -l`
        
        molRatio=2
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
            -nmol $numTBA \
            -o tert_box.gro >> $logFile 2>> $errFile 
        check tert_box.gro 

        gmx insert-molecules -ci tip3p.pdb \
            -f tert_box.gro \
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

        sed 's/POSRES/POSRES_IONS/' neutral_Ion3.itp > temp.itp 
        mv temp.itp neutral_Ion3.itp

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

production(){
    printf "\t\tProduction run............................" 
    if [ ! -f Production/$MOLEC.nopbc.gro ] ; then 
        create_dir Production
        
        cp Solvent_npt/neutral.top Production/.
        cp Solvent_npt/solvent_npt.gro Production/.
        cp Solvent_npt/*.itp Production/. 
        cd Production

        if [ ! -f $MOLEC.gro ] ; then 
            if [ ! -f $MOLEC.tpr ] ; then 
                gmx grompp -f $MDP/production_solvent.mdp \
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
                -s $MOLEC.tpr \
                -center \
                -ur rect \
                -pbc mol \
                -o $MOLEC.nopbc.xtc >> $logFile 2>> $errFile 
            fi 
        check $MOLEC.nopbc.xtc 

        if [ ! -f $MOLEC.nopbc.gro ] ; then 
            echo 'Protein System' | gmx trjconv -f $MOLEC.gro \
                -s $MOLEC.tpr \
                -center \
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
        echo "r TBUT & a O" > selection.dat 
        echo "r SOL & a OW" >> selection.dat 
        echo "a CD1 or CD2 & r LEU" >> selection.dat 
        echo "r LYS & a NZ" >> selection.dat 
        echo "q" >> selection.dat 

        cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.gro \
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
#dssp 
#rgyr 
#minimage
rdf 
cd ../

printf "\n\n\t\t*** Program Ending    ***\n\n" 
