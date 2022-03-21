#!/bin/bash

#PBS -N KDtree_MPI
#PBS -q dssc
#PBS -l nodes=1:ppn=24
#PBS -l walltime=01:00:00
#PBS -o KDtree_MPI_last2_run4.out
#PBS -j oe

cd ${PBS_O_WORKDIR}

module load openmpi-4.1.1+gnu-9.3.0


if [ -f KDtree_MPI_last2 ]
then make clean
fi

make CFLAGS="-Wall -O3 -march=native -std=c11 -DNPS=1000000000"
#make USER_CFLAGS="-DNPS=100000000"

narr=(1 2 8 16 24)
for i in ${narr[@]};
do
   for j in {0..4}
   do  
       mpirun -np $i --map-by ppr:1:core ./KDtree_MPI_last2 >> data/output_run4/${i}_${j}
   done
done
  


