#!/bin/bash
## Jeremy First

TOP=${PWD}
MDP=$TOP/mdp_files
FF=$TOP/GMXFF
forceField=oplsaa
dim=7.600
fileName=$TOP/StartingStructures/folded.pdb
totSimTime=50
verbose=false
analysis=false

usage(){
    echo "USAGE: $0 -f prep/folded/unfolded [ options ]" 
    exit 
}

HELP(){
    echo 
    echo "This program runs molecular dynamic simulations of a peptide in water "
    echo 
    echo "Usage: $0 [options] "
    echo "  -f   folded, unfolded, or prep. Mandatory" 
    echo "  -c   Starting structure  (folded) PDB file: Default = StartingStructures/folded.pdb"
    echo "  -t   Maximum simulation time (ns) : Default = 50 "
    echo "  -a   Perform analyis on trajectory : Default = no"
    echo "  -d   Box dimension (nm) : Default = $dim " 
    echo "  -m   Location of the mdp_files : Default = mdp_files"
    echo "  -p   Location of force field files. : Default = GMXFF"
    echo "  -n   Name of force field : Default = oplsaa"
    echo "  -v   Print all options and quit." 
    echo "  -h   print this usage and exit "
    echo "" 
    exit
}

while getopts :f:c:t:d:am:p:n:vh opt; do 
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
    checkInput
    printf "\n\t\t*** Program Beginning $SOL_$fold $totSimTime (ns)***\n\n" 
    if [ ! -d $SOL ] ; then mkdir $SOL ; fi 
    cd $SOL
    if [ $fold == "prep" ] ; then 
        prep
    else  
        #production 
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
    #rgyr
    #minimage
    #rdf
    nopbc
    rmsd
    #cd_spectra 
    #cluster
    good-turing
    hbond
    volume
    excluded_volume
    cd ../
}

checkInput(){
    SOL=water
    MOLEC=${fold}_${SOL}
    logFile=$TOP/$SOL/$fold/$fold.log
    errFile=$TOP/$SOL/$fold/$fold.err
    if $verbose ; then 
        echo "MOLEC : $MOLEC" 
        echo "SOL : $SOL" 
        echo "Folded state : $fold" 
        echo "Input file name: $fileName"
        echo "Max simultaiton time: $totSimTime"
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
        cp Solvate/neutral_*.itp Solvent_steep/. 
        cp Solvate/posre_*.itp Solvent_steep/. 
        cd Solvent_steep

        gmx grompp -f $MDP/solvent_steep_easy.mdp \
            -p neutral.top \
            -c neutral.gro \
            -o solvent_steep_easy.tpr >> $logFile 2>> $errFile 
        check solvent_steep_easy.tpr 

        gmx mdrun -deffnm solvent_steep_easy -nt 1 >> $logFile 2>> $errFile 
        check solvent_steep_easy.gro 
        
        gmx grompp -f $MDP/solvent_steep.mdp \
            -p neutral.top \
            -c solvent_steep_easy.gro \
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
                    ibrun mdrun_mpi -deffnm $MOLEC \
                        -s $simTime.tpr \
                        -cpi $MOLEC.cpt >> $logFile 2>> $errFile  
                else 
                    ibrun mdrun_mpi -deffnm $MOLEC \
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

        if [ ! -f scount.xvg ] ; then 
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
        fi 
        check scount.xvg ss.xpm area.xpm 

        if [ ! -f area.eps ] ; then 
            gmx xpm2ps -f area.xpm \
                -di $TOP/m2p_files/ps.m2p \
                -by 10 \
                -o area.eps >> $logFile 2>> $errFile 
        fi 
        check area.eps
        if [ ! -f area.png ] ; then 
            ps2pdf area.eps 
            convert area.pdf area.png 
        fi 

        if [ ! -f ss.eps ] ; then 
            gmx xpm2ps -f ss.xpm \
                -di $TOP/m2p_files/ps.m2p \
                -by 10 \
                -o ss.eps >> $logFile 2>> $errFile 
        fi 
        check ss.eps

        if [ ! -f ss.png ] ; then 
            ps2pdf ss.eps 
            convert ss.pdf ss.png 
        fi 

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
        if 'Turn' in line : 
            turn =line.split()[1][1:] 
        if 'B-Sheet' in line : 
            bsheet = line.split()[1][1:] 
        if 'B-Bridge' in line : 
            bbridge= line.split()[1][1:] 
        if 'Bend' in line : 
            bend = line.split()[1][1:] 
header -= 2

data = np.genfromtxt('scount.xvg',skip_header=header)
col1 = data[:,0]

try : 
    col2a = data[:,int(alpha)+1]
except NameError : 
    col2a = np.zeros(len(data[:,0])) 
try : 
    col2b = data[:,int(three)+1]
except NameError : 
    col2b = np.zeros(len(data[:,0])) 
try : 
    col2c = data[:,int(five)+1]
except NameError : 
    col2c = np.zeros(len(data[:,0])) 
try : 
    col2d = data[:,int(turn)+1]
except NameError : 
    col2d = np.zeros(len(data[:,0])) 
col2=col2a + col2b + col2c + col2d

try : 
    col3a = data[:,int(bsheet)+1]
except NameError : 
    col3a = np.zeros(len(data[:,0]))  
try : 
    col3b = data[:,int(bbridge)+1]
except NameError : 
    col3b = np.zeros(len(data[:,0]))  
try : 
    col3c = data[:,int(bend)+1]
except NameError : 
    col3c = np.zeros(len(data[:,0]))  
col3=col3a + col3b + col3c 

col4 = 20 - col2 - col3

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

        gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -ref 'name CD1 CD2 and resname LEU' \
            -sel 'resname SOL and name OW' \
            -o leu_wat.xvg >> $logFile 2>> $errFile 
        check leu_wat.xvg 

        gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -ref 'name NZ and resname LYS' \
            -sel 'resname SOL and name OW' \
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

        echo "Backbone Backbone" | gmx rms -f ../nopbc/nopbc.xtc \
            -s $fileName \
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

cluster(){
    printf "\t\tCluster analysis.........................." 
    if [ ! -d ../both ] ; then 
        cd ../
        create_dir both
        cd both 
        clean 
        cd ../$fold
    fi 

    cd ../both
    if [ ! -f trajcat/both.xtc ] ; then 
        create_dir trajcat 
        cd trajcat 
        clean 

        if [[ ! -f ../../folded/Production/folded_$SOL.xtc || ! -f ../../unfolded/Production/unfolded_$SOL.xtc ]] ; then 
            echo "Cannot create concatenated trajectory without production trajectories." 
            exit 
        fi 

        gmx trjcat -f ../../folded/Production/folded_$SOL.xtc ../../unfolded/Production/unfolded_$SOL.xtc \
            -cat \
            -b 625000 \
            -dt 100 \
            -o both.xtc >> $logFile 2>> $errFile
        check both.xtc 
        cd ../
    fi 

    if [ ! -f cluster/clust-size.xvg ] ; then 
        create_dir cluster
        cd cluster
        clean 

        echo "Backbone Backbone" | gmx rms -s ../../folded/Production/folded_$SOL.tpr \
            -f ../trajcat/both.xtc \
            -m rmsd.xpm >> $logFile 2>> $errFile 
        check rmsd.xpm 

        echo "Backbone Protein" | gmx cluster -s ../../folded/Production/folded_$SOL.tpr \
            -f ../trajcat/both.xtc \
            -dm rmsd.xpm \
            -cutoff 0.40 \
            -method gromos \
            -tr clust-trans.xpm \
            -cl clusters.pdb \
            -sz clust-size.xvg >> $logFile 2>> $errFile 
        check clust-size.xvg clusters.pdb 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
    cd ../$fold
}

good-turing(){
    printf "\t\tGood-Turing Stats........................." 
    gtScripts=$TOP/good-turing_scripts
    if [ ! -f good-turing/good_turing.rmsd.tar ] ; then 
        create_dir good-turing
        cd good-turing
        clean 

        sampling=1000 ##every 1 ns 
        equilTime=600000 ##600 ns 
        if [ ! -s rmsd.xpm ] ; then 
            echo "Backbone Backbone" | gmx rms -s ../Production/$MOLEC.tpr \
                -f ../Production/$MOLEC.xtc \
                -dt $sampling \
                -b $equilTime \
                -m rmsd.xpm >> $logFile 2>> $errFile 
        fi 
        check rmsd.xpm 

        python $gtScripts/xpm2dat.py rmsd.xpm 
        check rmsd.dat 

        ##https://github.com/pkoukos/GoodTuringMD
        ##  Adapted to avoid interacitve file name input
        echo "source('$gtScripts/Good_Turing.R')" > run_gt_stat.R

        RScript --no-save --no-restore --verbose run_gt_stat.R > good-turing.log 2>&1 
        check good_turing.rmsd.tar

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

hbond(){
    equilTime=600  ## ns 
    sampling=10    ## ns 
    printf "\t\tAnalyzing hbonds to peptide..............." 
    if [[ ! -f hbond/sidechain.xvg || ! -f hbond/nearby_sidechain.xvg \
        || ! -f hbond/mainchain.xvg || ! -f hbond_nearby_mainchain.xvg \
        || ! -f hbond/protein.xvg || ! -f hbond/nearby_protein.xvg ]] ; then 
        create_dir hbond
        cd hbond
        clean 

        if [ ! -f protein.xvg ] ; then 
            echo 'Protein Water' | gmx hbond -s ../Production/$MOLEC.tpr \
                -f ../Production/$MOLEC.xtc \
                -b $((equilTime*1000)) \
                -dt $((sampling*1000)) \
                -num protein.xvg >> $logFile 2>> $errFile 
        fi 
        check protein.xvg 

        if [ ! -f nearby_protein.xvg ] ; then 
            echo 'name "OW" and same residue as within 0.27 of group "Protein"' > selection.dat 
            cat selection.dat | gmx select -s ../Production/$MOLEC.tpr \
                -f ../Production/$MOLEC.xtc \
                -b $((equilTime*1000)) \
                -dt $((sampling*1000)) \
                -os nearby_protein.xvg >> $logFile 2>> $errFile 
        fi 
        check nearby_protein.xvg 

        if [ ! -f mainchain.xvg ] ; then 
            echo 'MainChain+H Water' | gmx hbond -s ../Production/$MOLEC.tpr \
                -f ../Production/$MOLEC.xtc \
                -b $((equilTime*1000)) \
                -dt $((sampling*1000)) \
                -num mainchain.xvg >> $logFile 2>> $errFile 
        fi 
        check mainchain.xvg 

        if [ ! -f nearby_mainchain.xvg ] ; then 
            echo 'name "OW" and same residue as within 0.27 of group "MainChain+H"' > selection.dat
            cat selection.dat | gmx select -s ../Production/$MOLEC.tpr \
                -f ../Production/$MOLEC.xtc \
                -b $((equilTime*1000)) \
                -dt $((sampling*1000)) \
                -os nearby_mainchain.xvg >> $logFile 2>> $errFile 
        fi 
        check nearby_mainchain.xvg 

        if [ ! -f sidechain.xvg ] ; then 
            echo 'SideChain Water' | gmx hbond -s ../Production/$MOLEC.tpr \
                -f ../Production/$MOLEC.xtc \
                -b $((equilTime*1000)) \
                -dt $((sampling*1000)) \
                -num sidechain.xvg >> $logFile 2>> $errFile 
        fi
        check sidechain.xvg 

        if [ ! -f nearby_sidechain.xvg ] ; then 
            echo 'name "OW" and same residue as within 0.27 of group "SideChain"' > selection.dat
            cat selection.dat | gmx select -s ../Production/$MOLEC.tpr \
                -f ../Production/$MOLEC.xtc \
                -b $((equilTime*1000)) \
                -dt $((sampling*1000)) \
                -os nearby_sidechain.xvg >> $logFile 2>> $errFile 
        fi 
        check nearby_sidechain.xvg 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

volume(){
    equilTime=600  ## ns 
    printf "\t\tCalculating volume of peptide............." 
    if [[ ! -f volume/volume.xvg || ! -f volume/vdw.xvg ]] ; then 
        create_dir volume
        cd volume
        clean 

        if [ ! -f volume.xvg ] ; then 
            echo 'Protein' | gmx sasa -s ../Production/$MOLEC.tpr \
                -f ../Production/$MOLEC.xtc \
                -b $((equilTime*1000)) \
                -ndots 240 \
                -tv volume.xvg >> $logFile 2>> $errFile 
        fi 
        check volume.xvg 

        if [ ! -f vdw.xvg ] ; then 
            echo 'Protein' | gmx sasa -s ../Production/$MOLEC.tpr \
                -f ../Production/$MOLEC.xtc \
                -b $((equilTime*1000)) \
                -ndots 240 \
                -probe 0 \
                -tv vdw.xvg >> $logFile 2>> $errFile 
        fi 
        check vdw.xvg 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

excluded_volume(){
    equilTime=600  ## ns 
    sampling=10    ## ns 
    printf "\t\tExcluded volume of peptide + water:      \n"

    create_dir excluded_volume
    cd excluded_volume
    clean 

    printf "\t\t     Pre-compressing trajectory..........."
    if [ ! -f whole.xtc ] ; then 
        echo 'Protein System' | gmx trjconv -s ../Production/$MOLEC.tpr \
            -f ../Production/$MOLEC.xtc \
            -pbc mol \
            -b $((equilTime*1000)) \
            -dt $((sampling*1000)) \
            -center \
            -o whole.xtc >> $logFile 2>> $errFile 
        check whole.xtc 
        printf "Success\n" 
    else 
        printf "Skipped\n"
    fi 

    for radius in `seq -w 0.05 0.02 0.50` ; do 
        printf "\t\t     Radius: %4s........................." $radius 

        if [ ! -f size_${radius}.xvg ] ; then 
            gmx select -s ../Production/$MOLEC.tpr \
                -f whole.xtc \
                -b $((equilTime*1000)) \
                -dt $((sampling*1000)) \
                -select "name OW and group \"Water\" and same residue as within $radius of group \"Protein\" " \
                -os size_${radius}.xvg >> $logFile 2>> $errFile 
        fi 
        check size_${radius}.xvg 

        for probe in `seq -w 0.00 0.01 0.14` ; do 
            if [ ! -f volume_p_${probe}_r_${radius}.xvg ] ; then 
                echo "group \"Protein\" or (group \"Water\" and same residue as within $radius of group \"Protein\")" > selection.dat 

                cat selection.dat |  gmx sasa \
                    -s ../Production/$MOLEC.tpr \
                    -f whole.xtc \
                    -b $((equilTime*1000)) \
                    -dt $((sampling*1000)) \
                    -probe ${probe} \
                    -ndots 240 \
                    -o area_p_${probe}_r_${radius}.xvg \
                    -tv volume_p_${probe}_r_${radius}.xvg >> $logFile 2>> $errFile
            fi 
            check volume_p_${probe}_r_${radius}.xvg
        done 

        printf "Success\n" 
    done 
    cd ../
}

main
