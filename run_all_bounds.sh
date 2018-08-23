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
#    echo "    -p All prep    " 
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
#       p) 
#         prep=true
#         ;; 
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

solvent="sam" 
molList="folded unfolded" 

if $folded ; then 
    for res in "0" "1" ; do
        sol="sam"

        printf "\n\n\t\t*** $solvent $res restraints ***\n\n" 

        if ! $stampede ; then 
            bash run_bounds_${solvent}.sh -f folded -n $res
        else 
            echo "#!/bin/bash" > submit_f_${sol}_${res}
            echo >> submit_f_${sol}_${res}
            echo "#SBATCH -J f_${sol}_${res} " >> submit_f_${sol}_${res}
            echo "#SBATCH -o f_${sol}_${res}.o%j" >> submit_f_${sol}_${res}
	    echo "#SBATCH -N 1" >> submit_f_${sol}_${res}
            echo "#SBATCH -n 48 " >> submit_f_${sol}_${res}
            echo "#SBATCH -p skx-normal " >> submit_f_${sol}_${res}
            echo "#SBATCH -t 48:00:00" >> submit_f_${sol}_${res}
            echo "#SBATCH -A Ras" >> submit_f_${sol}_${res}
            echo "#SBATCH --mail-user=jeremy_first@utexas.edu" >> submit_f_${sol}_${res}
            echo "#SBATCH --mail-type=all" >> submit_f_${sol}_${res}
            
            echo >> submit_f_${sol}_${res}
            echo "module load gromacs " >> submit_f_${sol}_${res} 
                                                                           
            echo >> submit_f_${sol}_${res}
            echo "bash run_bounds_${solvent}.sh -f folded -n $res -t $simTime " >> submit_f_${sol}_${res}

            sbatch submit_f_${sol}_${res} 
            fi 
        done 
    fi 

if $unfolded ; then 
    for res in "0" "1" ; do
        sol="sam"

        printf "\n\n\t\t*** $solvent $res restraints ***\n\n" 

        if ! $stampede ; then 
            bash run_bounds_${solvent}.sh -f unfolded -n $res
        else 
            echo "#!/bin/bash" > submit_u_${sol}_${res}
            echo >> submit_u_${sol}_${res}
            echo "#SBATCH -J u_${sol}_${res} " >> submit_u_${sol}_${res}
            echo "#SBATCH -o u_${sol}_${res}.o%j" >> submit_u_${sol}_${res}
	    echo "#SBATCH -N 1" >> submit_u_${sol}_${res}
            echo "#SBATCH -n 48 " >> submit_u_${sol}_${res}
            echo "#SBATCH -p skx-normal " >> submit_u_${sol}_${res}
            echo "#SBATCH -t 48:00:00" >> submit_u_${sol}_${res}
            echo "#SBATCH -A Ras" >> submit_u_${sol}_${res}
            echo "#SBATCH --mail-user=jeremy_first@utexas.edu" >> submit_u_${sol}_${res}
            echo "#SBATCH --mail-type=all" >> submit_u_${sol}_${res}
            
            echo >> submit_u_${sol}_${res}
            echo "module load gromacs " >> submit_u_${sol}_${res} 
                                                                           
            echo >> submit_u_${sol}_${res}
            echo "bash run_bounds_${solvent}.sh -f unfolded -n $res -t $simTime" >> submit_u_${sol}_${res}

            sbatch submit_u_${sol}_${res} 
            fi 
        done 
    fi 

