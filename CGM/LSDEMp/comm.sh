#!/bin/bash

#SBATCH --job-name="iso_linomp" 
#SBATCH --output=res.txt
#SBATCH --export=ALL
#SBATCH -p oro
#SBATCH -w compute-0-[7]
#SBATCH -n 24

ulimit -l unlimited

export OMP_NUM_THREADS=24

#mpiCC -fopenmp -I ~/eigen -o mainISO MainIso.cpp

/usr/lib64/openmpi/bin/mpirun -np 1 -x OMP_NUM_THREADS=24 mainISO
