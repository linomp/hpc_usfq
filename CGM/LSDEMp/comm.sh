#!/bin/bash

#SBATCH --job-name="iso_compress_1000clones"
#SBATCH --output="iso_compress_1000clones.%j"
#SBATCH --export=ALL
#SBATCH -p oro
#SBATCH -w compute-0-[1]
#SBATCH -n 24

ulimit -l unlimited

export OMP_NUM_THREADS=24

#mpiCC -fopenmp -I ~/eigen -o mainISO MainIso.cpp

mpirun -np 1 -x OMP_NUM_THREADS=24 mainISO
