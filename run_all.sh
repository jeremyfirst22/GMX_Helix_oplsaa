#!/bin/bash

solventList="water tert sam" 
molList="folded unfolded" 
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
    for molec in $molList ; do 
        if [[ "$molec" == "folded" ]] ; then 
            mol="fol" 
        elif [[ "$molec" == "unfolded" ]] ; then 
            mol="unf" 
        else 
            echo "ERROR: Unknown molecule $molec" ; continue 
            fi 
        if [ ! -f StartingStructures/$molec.pdb ] ; then 
            echo "ERROR: $molec.pdb not found!" 
            continue 
            fi 

        printf "\n\n\t\t*** $molec ***\n\n" 

        if [ ! -f submit_${mol}_${sol} ] ; then 
            echo "#!/bin/bash" > submit_${mol}_${sol}
            echo >> submit_${mol}_${sol}
            echo "#SBATCH -J ${mol}_${sol} " >> submit_${mol}_${sol}
            echo "#SBATCH -o ${mol}_${sol}.o%j" >> submit_${mol}_${sol}
            echo "#SBATCH -n 16 " >> submit_${mol}_${sol}
            echo "#SBATCH -p normal " >> submit_${mol}_${sol}
            echo "#SBATCH -t 48:00:00" >> submit_${mol}_${sol}
            echo "#SBATCH -A Ras" >> submit_${mol}_${sol}
            echo "#SBATCH --mail-user=jeremy_first@utexas.edu" >> submit_${mol}_${sol}
            echo "#SBATCH --mail-type=all" >> submit_${mol}_${sol}
            
            echo >> submit_${mol}_${sol}
            echo "module load boost " >> submit_${mol}_${sol}
            echo "module load cxx11 " >> submit_${mol}_${sol}
            echo "module load gromacs " >> submit_${mol}_${sol} 
                                                                           
            echo >> submit_${mol}_${sol}
            echo "bash run_${solvent}.sh StartingStructures/$molec.pdb 100" >> submit_${mol}_${sol}
            fi 
        #sbatch submit_${mol}_${sol} 
        bash run_${solvent}.sh StartingStructures/$molec.pdb 

    
        done 
     done 

