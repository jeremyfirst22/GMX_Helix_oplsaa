#!/bin/bash

all=false
prep=false
folded=false
unfolded=false
stampede=false
simTime=50

HELP(){
    echo 
    echo "This programs is a wrapper for the run_sam.sh, run_water.sh, and run_tert.sh scripts." 
    echo 
    echo "Usage: $0 [options] " 
    echo "    -s Run on stampede using sbatch " 
#    echo "    -a All systems " 
    echo "    -p All prep    " 
    echo "    -f All three folded " 
    echo "    -u All three unfolded system " 
    echo "    -t Simulation time (Default = 50ns) " 
    echo "    -h Print this help message and exit " 
    echo 
    exit
}

while getopts sapfut:h opt ; do 
    case $opt in 
       s) 
         stampede=true 
         ;; 
#       a) 
#         all=true
#         ;; 
       p) 
         prep=true
         ;; 
       f)
         folded=true
         ;; 
       u) 
         unfolded=true
         ;; 
       t) 
         simTime=$OPTARG
         ;; 
       :) 
         echo " option $OPTARG requires an argument"
         usage
         ;; 
       h) 
         HELP 
         ;; 
       \?) 
         echo "Inavalid option -$OPTARG " 
         HELP 
         ;; 
     esac 
     done 

if ! $prep && ! $folded && ! $unfolded ; then 
    echo "Nothing to do! Exitting. " 
    echo "   Use $0 -h for more information" 
    exit 
    fi 

if [ ! -z $WORK ] && ! $stampede ; then 
    echo "Must use -s flag to run on stampede"
    exit
    fi 
if [ $simTime -lt 50 ] ; then 
    echo "Simulation time must be greater than 50!" 
    exit 
    fi 

solventList="water tert sam" 
molList="folded unfolded" 

if $prep ; then 
    for solvent in $solventList ; do 
        if [[ "$solvent" = "water" ]] ; then 
            sol="wat"
        elif [[ "$solvent" = "tert" ]] ; then 
            sol="tba"
        elif [[ "$solvent" = "sam" ]] ; then 
            sol="sam"
        else 
            echo "ERROR: Unknown solvent $solvent" ; continue 
            fi 

        printf "\n\n\t\t*** $solvent ***\n\n" 

        if ! $stampede ; then 
            bash run_${solvent}.sh -f prep 
        else 
            echo "#!/bin/bash" > submit_prep_${sol}
            echo >> submit_prep_${sol}
            echo "#SBATCH -J prep_${sol} " >> submit_prep_${sol}
            echo "#SBATCH -o prep_${sol}.o%j" >> submit_prep_${sol}
            echo "#SBATCH -n 16 " >> submit_prep_${sol}
            echo "#SBATCH -p normal " >> submit_prep_${sol}
            echo "#SBATCH -t 03:00:00" >> submit_prep_${sol}
            echo "#SBATCH -A Ras" >> submit_prep_${sol}
            echo "#SBATCH --mail-user=jeremy_first@utexas.edu" >> submit_prep_${sol}
            echo "#SBATCH --mail-type=all" >> submit_prep_${sol}
            
            echo >> submit_prep_${sol}
            echo "module load boost " >> submit_prep_${sol}
            echo "module load cxx11 " >> submit_prep_${sol}
            echo "module load gromacs " >> submit_prep_${sol} 
                                                                           
            echo >> submit_prep_${sol}
            echo "bash run_${solvent}.sh -f prep " >> submit_prep_${sol}

            sbatch submit_prep_${sol} 
            fi 
        done 
    fi 

if $folded ; then 
    for solvent in $solventList ; do 
        if [[ "$solvent" = "water" ]] ; then 
            sol="wat"
        elif [[ "$solvent" = "tert" ]] ; then 
            sol="tba"
        elif [[ "$solvent" = "sam" ]] ; then 
            sol="sam"
        else 
            echo "ERROR: Unknown solvent $solvent" ; continue 
            fi 

        printf "\n\n\t\t*** $solvent ***\n\n" 

        if ! $stampede ; then 
            bash run_${solvent}.sh -f folded
        else 
            echo "#!/bin/bash" > submit_fol_${sol}
            echo >> submit_fol_${sol}
            echo "#SBATCH -J fol_${sol} " >> submit_fol_${sol}
            echo "#SBATCH -o fol_${sol}.o%j" >> submit_fol_${sol}
            echo "#SBATCH -n 16 " >> submit_fol_${sol}
            echo "#SBATCH -p normal " >> submit_fol_${sol}
            echo "#SBATCH -t 48:00:00" >> submit_fol_${sol}
            echo "#SBATCH -A Ras" >> submit_fol_${sol}
            echo "#SBATCH --mail-user=jeremy_first@utexas.edu" >> submit_fol_${sol}
            echo "#SBATCH --mail-type=all" >> submit_fol_${sol}
            
            echo >> submit_fol_${sol}
            echo "module load boost " >> submit_fol_${sol}
            echo "module load cxx11 " >> submit_fol_${sol}
            echo "module load gromacs " >> submit_fol_${sol} 
                                                                           
            echo >> submit_fol_${sol}
            echo "bash run_${solvent}.sh -f folded -t $simTime " >> submit_fol_${sol}

            sbatch submit_fol_${sol} 
            fi 
        done 
    fi 

if $unfolded ; then 
    for solvent in $solventList ; do 
        if [[ "$solvent" = "water" ]] ; then 
            sol="wat"
        elif [[ "$solvent" = "tert" ]] ; then 
            sol="tba"
        elif [[ "$solvent" = "sam" ]] ; then 
            sol="sam"
        else 
            echo "ERROR: Unknown solvent $solvent" ; continue 
            fi 

        printf "\n\n\t\t*** $solvent ***\n\n" 

        if ! $stampede ; then 
            bash run_${solvent}.sh -f unfolded
        else 
            echo "#!/bin/bash" > submit_unf_${sol}
            echo >> submit_unf_${sol}
            echo "#SBATCH -J unf_${sol} " >> submit_unf_${sol}
            echo "#SBATCH -o unf_${sol}.o%j" >> submit_unf_${sol}
            echo "#SBATCH -n 16 " >> submit_unf_${sol}
            echo "#SBATCH -p normal " >> submit_unf_${sol}
            echo "#SBATCH -t 48:00:00" >> submit_unf_${sol}
            echo "#SBATCH -A Ras" >> submit_unf_${sol}
            echo "#SBATCH --mail-user=jeremy_first@utexas.edu" >> submit_unf_${sol}
            echo "#SBATCH --mail-type=all" >> submit_unf_${sol}
            
            echo >> submit_unf_${sol}
            echo "module load boost " >> submit_unf_${sol}
            echo "module load cxx11 " >> submit_unf_${sol}
            echo "module load gromacs " >> submit_unf_${sol} 
                                                                           
            echo >> submit_unf_${sol}
            echo "bash run_${solvent}.sh -f unfolded -t $simTime" >> submit_unf_${sol}

            sbatch submit_unf_${sol} 
            fi 
        done 
    fi 

