#!/bin/bash
## Jeremy First

TOP=${PWD}
MDP=$TOP/mdp_files
FF=$TOP/GMXFF
SAM=$TOP/sam
forceField=oplsaa
saltConcentration=0
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
zBuff=0.5 #nm. This is the spacing of the SAM layer off the z intercept. Must not be zero, 
          ##    or the SAM splits during NPT box growth. 

verbose=false
analysis=false

usage(){
    echo "USAGE: $0 -f prep/folded/unfolded -n 0/1 [ options ]" 
    exit 
}

HELP(){
    echo 
    echo "This program runs molecular dynamic simulations of a peptide on a SAM surface" 
    echo 
    echo "Usage: $0 [options] "
    echo "  -f   folded, unfolded, or prep. Mandatory" 
    echo "  -n   Number of bounds : 1 or 0. Mandatory"
    echo "  -c   Starting structure  (folded) PDB file: Default = StartingStructures/folded.pdb"
    echo "  -t   Maximum simulation time (ns) : Default = 50 "
    echo "  -a   Perform analyis on trajectory : Default = no"
    echo "  -s   Salt conentration (mM) : Default = 0 (neutralizing)" 
    echo "  -n   Number of decanethiol molecules in each direction of the layer : Default = 14" 
    echo "  -d   Decanethiol structure PDB file : Default = StartingStructures/decanethiol.pdb" 
    echo "  -D   Distance (nm) between each decanethiol molecule : Default = 0.497" 
    echo "  -g   Distnace (nm) between the alpha carbon of glycines to C10 of the SAM layer" 
    echo "              Default = 0.6822 from Geoemetry optimization, Spartan '16, DFT wB97X-D, 6-311+g**"
    echo "  -r   Restraint force constant (kJ/mol/nm) of the "bond" from glycine to nearest decanethiol C10"
    echo "              Default = 1000 "
    echo "  -R   Restraint force constant (kJ/mol/nm) of the "bond" keeping sulfur atoms in plane" 
    echo "              Default = 200000" 
    echo "  -z   Distance of the bottom of the SAM layer to the z-axis of the simulation box. " 
    echo "              Default = 0.5 " 
    echo "  -m   Location of the mdp_files : Default = mdp_files"
    echo "  -p   Location of force field files. : Default = GMXFF"
    echo "  -n   Name of force field : Default = oplsaa"
    echo "  -v   Print all options and quit." 
    echo "  -h   print this usage and exit "
    echo "" 
    exit
}

while getopts :f:n:c:t:as:n:d:D:g:r:R:z:m:p:n:vh opt; do 
   case $opt in 
      f) 
        fold=$OPTARG
        ;; 
      n)
        numRestraints=$OPTARG
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
      s) 
        saltConcentration=$OPTARG
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
      z) 
        zBuff=$OPTARG
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

    if [ $numRestraints -eq 0 ] ; then
        SOL="not_bound_"$SOL
    elif [ $numRestraints -eq 1 ] ; then
        SOL="single_bound_"$SOL
    else
        echo "ERROR: Only not bound or single bound supported"
        exit
    fi
    MOLEC=folded_$SOL
    logFile=$TOP/$SOL/$fold/$fold.log
    errFile=$TOP/$SOL/$fold/$fold.err

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
    copy_from_sam_prep
    cd ../
}

analysis(){
    if [ ! -d $fold ] ; then mkdir $fold ; fi 
    cd $fold 
    dssp
    rgyr
    #minimage
    #rdf
    nopbc
    rmsd 
    cd_spectra
    #cluster
    good-turing
    order
    cd ../
}

checkInput(){
    if [ $saltConcentration -eq 0 ] ; then 
        SOL=sam
    else 
        SOL=sam_$saltConcentration
    fi 
    if $verbose ; then 
        echo "MOLEC : $MOLEC"
        echo "SOL : $SOL"
        echo "Folded state : $fold" 
        echo "Number of restraints : $numRestraints"
        echo "Input file name: $fileName"
        echo "Max simultaiton time: $totSimTime"
        echo "Salt concentration : $saltConcentration"
        echo "Number of decanethiols: $numMols" 
        echo "Decanethiol starting structure : $decanethiolFile" 
        echo "Distance between each decanethiol : $spacing"
        echo "Distance of glycine to SAM : $glyDist" 
        echo "Glycine restraint : $glyRest"
        echo "Sulfur restraint : $sulRest" 
        echo "Perform analysis : $analysis " 
        echo "Verbose = $verbose"
        echo "Distance to z-axis: $zBuff" 
        echo "mdp files = $MDP " 
        echo "sam run = $SAM "
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
    if [ ! -z $numRestraints ] ; then
        if [[ $numRestraints != '0' && $numRestraints != '1' ]] ; then
            echo "ERROR: $numRestraints not recognized"
            HELP
            fi
    else
        echo "ERROR: Must specify number of restraints 0 or 1"
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
    if [[ $saltConcentration < 0 ]] ; then 
        echo "ERROR: Negative salt conentration given" 
        exit 
        fi 
    if [[ ! -d $SAM && $fold == 'prep' ]] ; then
        echo ; echo "ERROR: sam run directory not found"
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

copy_from_sam_prep(){
    printf "\t\tPrepping from normal peptide on SAM......."
    ##Must copy from sam prep, since heating and cooling done with harmonic restraints.
    ## If heating, cooling done without harmonic restraints, unfolded simulation begins
    ## with peptide fully stretched away from SAM layer. Would take very long to simulate
    ## back to SAM layer, and has artifact of interacting with bottom of layer across PBC.
    if [[ ! -s Cooling/cooling.gro || ! -s System_npt/system_npt.gro || ! -s System_npt/system.top ]] ; then
        if [[ ! -s $SAM/prep/Cooling/cooling.gro || ! -s $SAM/prep/System_npt/system_npt.gro ]] ; then
            echo "ERROR: Must do run_sam.sh -f prep first! This program relies on normal prep of peptide on SAM"
            exit
        fi
        create_dir Cooling
        cd Cooling
        clean
        cp $SAM/prep/Cooling/cooling.gro .
        cd ../

        create_dir System_npt
        cd System_npt
        clean
        cp $SAM/prep/System_npt/system_npt.gro .


        if [ $numRestraints -eq 0 ] ; then
            #cp $SAM/prep/Build_system/neutral.top system.top ##before restraints added
            cp $SAM/prep/Build_system/system.top system.top ##with both restraints, must remove both
            cp $SAM/prep/Build_system/*.itp .

            cp distance_restraints.itp temp.itp
            sed '$ d' temp.itp | sed '$ d' > distance_restraints.itp
        elif [ $numRestraints -eq 1 ] ; then
            cp $SAM/prep/Build_system/system.top system.top ##with both restraints, must remove one
            cp $SAM/prep/Build_system/*.itp .

            cp distance_restraints.itp temp.itp
            sed '$ d' temp.itp > distance_restraints.itp
        else
            echo "ERROR: Only 1 or 0 restraints supported"
            exit
        fi

        clean
        printf "Success\n"
        cd ../
    else
        printf "Skipped\n"
        fi
}

production(){
    printf "\t\tProduction run............................" 
    if [[ ! -s prep/Cooling/cooling.gro || ! -s prep/System_npt/system_npt.gro ]] ; then 
        echo ; echo 
        echo "ERROR: You must run prep first!" 
        fi 
    check prep/Cooling/cooling.gro prep/System_npt/system_npt.gro 

    if [ ! -d $fold ] ; then mkdir $fold ; fi 
    cd $fold 
    MOLEC=${fold}_${SOL}

    if [ ! -f Production/${MOLEC}_${totSimTime}ns.gro ] ; then 
        printf "\n" 
        create_dir Production
        
        cp ../prep/System_npt/system.top Production/.
        cp ../prep/System_npt/*.itp Production/. 

        if [ $fold = "folded" ] ; then 
            cp ../prep/System_npt/system_npt.gro Production/.
            startStructure=system_npt.gro 
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
                    ibrun mdrun_mpi -deffnm $MOLEC \
                        -s $simTime.tpr \
                        -cpi $MOLEC.cpt >> $logFile 2>> $errFile  
                else 
                    #ibrun mdrun_mpi -deffnm $MOLEC \
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
    if [ ! -f rdf/sam_lys.xvg ] ; then 
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

        gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -ref 'name CD1 CD2 and resname LEU' \
            -sel 'name C10 and resname LIG' \
            -o sam_leu.xvg >> $logFile 2>> $errFile 
        check sam_leu.xvg 

        gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -ref 'name NZ and resname LYS' \
            -sel 'name C10 and resname LIG' \
            -o sam_lys.xvg >> $logFile 2>> $errFile 
        check sam_lys.xvg 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

nopbc(){
    printf "\t\tCreating nopbc copy......................." 
    if [ ! -f nopbc/fit.xtc ] ; then 
        create_dir nopbc
        cd nopbc
        clean 

        if [ $fold = "folded" ] ; then 
            startStructure=system_npt.gro 
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

        echo "1 | 14" > selection.dat 
        echo "q" >> selection.dat 

        cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.tpr >> $logFile 2>> $errFile 
        check index.ndx 

        echo "Protein_LIG" | gmx trjconv -f ../Production/$startStructure \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -pbc mol \
            -ur compact \
            -o fit.gro >> $logFile 2>> $errFile 
        check fit.gro 

        echo "Protein Protein_LIG" | gmx trjconv -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -fit transxy \
            -pbc res \
            -dt 40 \
            -ur compact \
            -o fit.xtc >> $logFile 2>> $errFile 
        check fit.xtc 

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
    if [ ! -f cluster/clust-size.xvg ] ; then 
        create_dir cluster
        cd cluster
        clean 

        echo "Backbone Backbone" | gmx rms -s ../Production/$MOLEC.tpr \
            -f ../Production/$MOLEC.xtc \
            -dt 100 \
            -m rmsd.xpm >> $logFile 2>> $errFile 
        check rmsd.xpm 

        echo "Backbone Protein" | gmx cluster -s ../Production/$MOLEC.tpr \
            -f ../Production/$MOLEC.xtc \
            -dm rmsd.xpm \
            -cutoff 0.25 \
            -method gromos \
            -dt 100 \
            -tr clust-trans.xpm \
            -cl clusters.pdb \
            -sz clust-size.xvg >> $logFile 2>> $errFile 
        check clust-size.xvg clusters.pdb 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

good-turing(){
    printf "\t\tGood-Turing Stats........................." 
    gtScripts=$TOP/good-turing_scripts
    if [ ! -f good-turing/good_turing.rmsd.tar ] ; then 
        create_dir good-turing
        cd good-turing
        clean 

        sampling=1000 ##every 1 ns 
        equilTime=50000 ## 50 ns 
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

order(){
    printf "\t\tCalculating order of SAM.................." 
    if [ ! -f order/density.xvg ] ; then 
        create_dir order
        cd order
        clean 

        echo "r LIG & a S1 C1 C2 C3 C4 C5 C6 C7 C8 C9 C10" > selection.dat 
        echo "q" >> selection.dat 

        touch empty.ndx 
        cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.tpr \
            -n empty.ndx \
            -o heavyLIG.ndx >> $logFile 2>> $errFile 
        check heavyLIG.ndx 

        gmx density -s ../Production/$MOLEC.tpr \
            -f ../Production/$MOLEC.xtc \
            -n heavyLIG.ndx \
            -d Z \
            -sl 1000 \
            -e 1000 \
            -o density.xvg >> $logFile 2>> $errFile 
        check density.xvg

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

main
