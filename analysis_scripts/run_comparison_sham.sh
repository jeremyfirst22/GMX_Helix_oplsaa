#!/bin/bash

TOP=${PWD}

check(){
   for arg in $@ ; do
        if [ ! -s $arg ] ; then
            echo ; echo "ERROR: $arg missing. Exitting"
            exit
            fi
        done
}


if [ ! -d sham_comparison ] ; then mkdir sham_comparison ; fi 
cd sham_comparison 

if [ ! -d ordered ] ; then mkdir ordered ; fi 
cd ordered 

logFile=${TOP}/sham_comparison/ordered/sham.log 
errFile=${TOP}/sham_comparison/ordered/sham.err 

equilTime=150 

cat ../../single_bound_sam/folded/rgyr/gyrate.xvg > data.xvg 
cat ../../single_bound_sam/unfolded/rgyr/gyrate.xvg >> data.xvg 
cat ../../not_bound_sam/folded/rgyr/gyrate.xvg >> data.xvg 
cat ../../not_bound_sam/unfolded/rgyr/gyrate.xvg >> data.xvg 
echo '&' >> data.xvg 
cat ../../single_bound_sam/folded/rmsd/rmsd.xvg >> data.xvg 
cat ../../single_bound_sam/unfolded/rmsd/rmsd.xvg >> data.xvg 
cat ../../not_bound_sam/folded/rmsd/rmsd.xvg >> data.xvg 
cat ../../not_bound_sam/unfolded/rmsd/rmsd.xvg >> data.xvg 

gmx sham -f data.xvg \
    -ls gibbs.xpm \
    -lp prob.xpm \
    -n 2 \
    -b $((equilTime*1000)) \
    -xmin 0.75 0.150 0 \
    -xmax 1.75 1.100 0 \
    -nlevels 25  \
    -ngrid 100 >> $logFile 2>> $errFile 
check gibbs.xpm 

if [ -f gibbs.dat ] ; then 
    rm gibbs.dat 
fi 

python ../../analysis_scripts/xpm2dat.py gibbs.xpm 
check gibbs.dat 

cd ../

if [ ! -d disordered ] ; then mkdir disordered ; fi 
cd disordered 

logFile=${TOP}/sham_comparison/disordered/sham.log 
errFile=${TOP}/sham_comparison/disordered/sham.err 

equilTime=150 

cat ../../disordered_single_bound_sam/folded/rgyr/gyrate.xvg > data.xvg 
cat ../../disordered_single_bound_sam/unfolded/rgyr/gyrate.xvg >> data.xvg 
cat ../../disordered_not_bound_sam/folded/rgyr/gyrate.xvg >> data.xvg 
cat ../../disordered_not_bound_sam/unfolded/rgyr/gyrate.xvg >> data.xvg 
echo '&' >> data.xvg 
cat ../../disordered_single_bound_sam/folded/rmsd/rmsd.xvg >> data.xvg 
cat ../../disordered_single_bound_sam/unfolded/rmsd/rmsd.xvg >> data.xvg 
cat ../../disordered_not_bound_sam/folded/rmsd/rmsd.xvg >> data.xvg 
cat ../../disordered_not_bound_sam/unfolded/rmsd/rmsd.xvg >> data.xvg 

gmx sham -f data.xvg \
    -ls gibbs.xpm \
    -lp prob.xpm \
    -n 2 \
    -b $((equilTime*1000)) \
    -xmin 0.75 0.150 0 \
    -xmax 1.75 1.100 0 \
    -nlevels 25  \
    -ngrid 100 >> $logFile 2>> $errFile 
check gibbs.xpm 

if [ -f gibbs.dat ] ; then 
    rm gibbs.dat 
fi 

python ../../analysis_scripts/xpm2dat.py gibbs.xpm 
check gibbs.dat 

cd ../

cd ../
        


    
