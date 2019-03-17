#!/bin/bash

#SBATCH --job-name="iso_linomp" 
#SBATCH --output=res.txt
#SBATCH --export=ALL
#SBATCH -p oro
#SBATCH -w compute-0-[1]
#SBATCH -n 16

ulimit -l unlimited

export OMP_NUM_THREADS=16

#mpiCC -fopenmp -I ~/eigen -o mainISO MainIso.cpp

mpirun -np 1 -x OMP_NUM_THREADS=16 mainISO
