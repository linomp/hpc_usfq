#!/bin/bash
#
#SBATCH --job-name=cuda_linomp
#SBATCH --output=res.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 16
#SBATCH --nodelist compute-0-12
#SBATCH --time=1:00


srun gcc -std=c99 -fopenmp evaluar_func.c -o ahmdal.exe -lm

export OMP_NUM_THREADS=2
srun ./ahmdal.exe 2

export OMP_NUM_THREADS=4
srun ./ahmdal.exe 4

export OMP_NUM_THREADS=8
srun ./ahmdal.exe 8

export OMP_NUM_THREADS=16
srun ./ahmdal.exe 16

export OMP_NUM_THREADS=32
srun ./ahmdal.exe 32
<<<<<<< HEAD
=======

>>>>>>> 865e6c0a1f920e85e9294e2b79445f6d741ea5c2
