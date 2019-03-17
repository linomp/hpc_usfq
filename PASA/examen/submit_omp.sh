#!/bin/bash
#
#SBATCH --job-name=cuda_linomp
#SBATCH --output=res.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 16
#SBATCH --nodelist compute-0-12
#SBATCH --time=1:00


srun gcc -std=c99 -fopenmp suma_vec.c -o suma.exe -lm

export OMP_NUM_THREADS=2
srun ./suma.exe 2

export OMP_NUM_THREADS=4
srun ./suma.exe 4

export OMP_NUM_THREADS=8
srun ./suma.exe 8

export OMP_NUM_THREADS=16
srun ./suma.exe 16 
