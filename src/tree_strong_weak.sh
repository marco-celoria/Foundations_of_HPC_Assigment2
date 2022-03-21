#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l walltime=1:00:00
#PBS -q dssc_gpu

cd /u/dssc/mceloria/kdtree2
module load openmpi-4.1.1+gnu-9.3.0
make code=serial precision=double
make code=openmp precision=double
make code=mpi precision=double

dir=/u/dssc/mceloria/kdtree3/data
if [[ ! -e $dir ]]; then
    mkdir $dir
elif [[ ! -d $dir ]]; then
    echo "$dir already exists but is not a directory" 1>&2
fi


echo "Strong scalability"                1> "${dir}/strong_scalability.dat"
echo "::::::::::::::::::::::::::::"      1>> "${dir}/strong_scalability.dat"

echo "Serial"                            1>> "${dir}/strong_scalability.dat"
echo "N=100000000"                       1>> "${dir}/strong_scalability.dat"
./tree_serial.x 100000000                1>> "${dir}/strong_scalability.dat"

echo "::::::::::::::::::::::::::::"      1>> "${dir}/strong_scalability.dat"
echo "OpenMP"                            1>> "${dir}/strong_scalability.dat"

export OMP_NUM_THREADS=1
echo "OMP_NUM_THREADS=1"                 1>> "${dir}/strong_scalability.dat"
echo "N=100000000"                       1>> "${dir}/strong_scalability.dat"
./tree_openmp.x 100000000                1>> "${dir}/strong_scalability.dat"
echo "============================"      1>> "${dir}/strong_scalability.dat"

export OMP_NUM_THREADS=2
echo "OMP_NUM_THREADS=2"                 1>> "${dir}/strong_scalability.dat"
echo "N=100000000"                       1>> "${dir}/strong_scalability.dat"
./tree_openmp.x 100000000                1>> "${dir}/strong_scalability.dat"
echo "============================"      1>> "${dir}/strong_scalability.dat"

export OMP_NUM_THREADS=4
echo "OMP_NUM_THREADS=4"                 1>> "${dir}/strong_scalability.dat"
echo "N=100000000"                       1>> "${dir}/strong_scalability.dat"
./tree_openmp.x 100000000                1>> "${dir}/strong_scalability.dat"
echo "============================"      1>> "${dir}/strong_scalability.dat"

export OMP_NUM_THREADS=8
echo "OMP_NUM_THREADS=8"                 1>> "${dir}/strong_scalability.dat"
echo "N=100000000"                       1>> "${dir}/strong_scalability.dat"
./tree_openmp.x 100000000                1>> "${dir}/strong_scalability.dat"
echo "============================"      1>> "${dir}/strong_scalability.dat"

export OMP_NUM_THREADS=16
echo "OMP_NUM_THREADS=16"                1>> "${dir}/strong_scalability.dat"
echo "N=100000000"                       1>> "${dir}/strong_scalability.dat"
./tree_openmp.x 100000000                1>> "${dir}/strong_scalability.dat"
echo "============================"      1>> "${dir}/strong_scalability.dat"

echo "::::::::::::::::::::::::::::"      1>> "${dir}/strong_scalability.dat"
echo "MPI"                               1>> "${dir}/strong_scalability.dat"

echo "MPI -np=1"                         1>> "${dir}/strong_scalability.dat"
echo "N=100000000"                       1>> "${dir}/strong_scalability.dat"
mpirun -np 1 ./tree_mpi.x 100000000      1>> "${dir}/strong_scalability.dat"
echo "============================"      1>> "${dir}/strong_scalability.dat"

echo "MPI -np=2"                         1>> "${dir}/strong_scalability.dat"
echo "N=100000000"                       1>> "${dir}/strong_scalability.dat"
mpirun -np 2 ./tree_mpi.x 100000000      1>> "${dir}/strong_scalability.dat"
echo "============================"      1>> "${dir}/strong_scalability.dat"

echo "MPI -np=4"                         1>> "${dir}/strong_scalability.dat"
echo "N=100000000"                       1>> "${dir}/strong_scalability.dat"
mpirun -np 4 ./tree_mpi.x 100000000      1>> "${dir}/strong_scalability.dat"
echo "============================"      1>> "${dir}/strong_scalability.dat"

echo "MPI -np=8"                         1>> "${dir}/strong_scalability.dat"
echo "N=100000000"                       1>> "${dir}/strong_scalability.dat"
mpirun -np 8 ./tree_mpi.x 100000000      1>> "${dir}/strong_scalability.dat"
echo "============================"      1>> "${dir}/strong_scalability.dat"

echo "MPI -np=16"                        1>> "${dir}/strong_scalability.dat"
echo "N=100000000"                       1>> "${dir}/strong_scalability.dat"
mpirun -np 16 ./tree_mpi.x 100000000     1>> "${dir}/strong_scalability.dat"
echo "============================"      1>> "${dir}/strong_scalability.dat"


echo "Weak scalability"                  1> "${dir}/weak_scalability.dat"
echo "::::::::::::::::::::::::::::"      1>> "${dir}/weak_scalability.dat"

echo "Serial"                            1>> "${dir}/weak_scalability.dat"
echo "N=100000000"                       1>> "${dir}/weak_scalability.dat"
./tree_serial.x 100000000                1>> "${dir}/weak_scalability.dat"

echo "::::::::::::::::::::::::::::"      1>> "${dir}/weak_scalability.dat"
echo "OpenMP"                            1>> "${dir}/weak_scalability.dat"

export OMP_NUM_THREADS=1
echo "OMP_NUM_THREADS=1"                 1>> "${dir}/weak_scalability.dat"
echo "N=100000000"                       1>> "${dir}/weak_scalability.dat"
./tree_openmp.x 100000000                1>> "${dir}/weak_scalability.dat"
echo "============================"      1>> "${dir}/weak_scalability.dat"

export OMP_NUM_THREADS=2
echo "OMP_NUM_THREADS=2"                 1>> "${dir}/weak_scalability.dat"
echo "N=200000000"                       1>> "${dir}/weak_scalability.dat"
./tree_openmp.x 200000000                1>> "${dir}/weak_scalability.dat"
echo "============================"      1>> "${dir}/weak_scalability.dat"

export OMP_NUM_THREADS=4
echo "OMP_NUM_THREADS=4"                 1>> "${dir}/weak_scalability.dat"
echo "N=400000000"                       1>> "${dir}/weak_scalability.dat"
./tree_openmp.x 400000000                1>> "${dir}/weak_scalability.dat"
echo "============================"      1>> "${dir}/weak_scalability.dat"

export OMP_NUM_THREADS=8
echo "OMP_NUM_THREADS=8"                 1>> "${dir}/weak_scalability.dat"
echo "N=800000000"                       1>> "${dir}/weak_scalability.dat"
./tree_openmp.x 800000000                1>> "${dir}/weak_scalability.dat"
echo "============================"      1>> "${dir}/weak_scalability.dat"

export OMP_NUM_THREADS=16
echo "OMP_NUM_THREADS=16"                1>> "${dir}/weak_scalability.dat"
echo "N=1600000000"                      1>> "${dir}/weak_scalability.dat"
./tree_openmp.x 1600000000               1>> "${dir}/weak_scalability.dat"
echo "============================"      1>> "${dir}/weak_scalability.dat"

echo "::::::::::::::::::::::::::::"      1>> "${dir}/weak_scalability.dat"
echo "MPI"                               1>> "${dir}/weak_scalability.dat"

echo "MPI -np=1"                         1>> "${dir}/weak_scalability.dat"
echo "N=100000000"                       1>> "${dir}/weak_scalability.dat"
mpirun -np 1 ./tree_mpi.x 100000000      1>> "${dir}/weak_scalability.dat"
echo "============================"      1>> "${dir}/weak_scalability.dat"

echo "MPI -np=2"                         1>> "${dir}/weak_scalability.dat"
echo "N=200000000"                       1>> "${dir}/weak_scalability.dat"
mpirun -np 2 ./tree_mpi.x 200000000      1>> "${dir}/weak_scalability.dat"
echo "============================"      1>> "${dir}/weak_scalability.dat"

echo "MPI -np=4"                         1>> "${dir}/weak_scalability.dat"
echo "N=400000000"                       1>> "${dir}/weak_scalability.dat"
mpirun -np 4 ./tree_mpi.x 400000000      1>> "${dir}/weak_scalability.dat"
echo "============================"      1>> "${dir}/weak_scalability.dat"

echo "MPI -np=8"                         1>> "${dir}/weak_scalability.dat"
echo "N=800000000"                       1>> "${dir}/weak_scalability.dat"
mpirun -np 8 ./tree_mpi.x 800000000      1>> "${dir}/weak_scalability.dat"
echo "============================"      1>> "${dir}/weak_scalability.dat"

echo "MPI -np=16"                        1>> "${dir}/weak_scalability.dat"
echo "N=1600000000"                      1>> "${dir}/weak_scalability.dat"
mpirun -np 16 ./tree_mpi.x 1600000000    1>> "${dir}/weak_scalability.dat"
echo "============================"      1>> "${dir}/weak_scalability.dat"

exit
    

