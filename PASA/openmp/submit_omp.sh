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
srun ./ahmdal.exe

export OMP_NUM_THREADS=4
srun ./ahmdal.exe

export OMP_NUM_THREADS=8
srun ./ahmdal.exe

export OMP_NUM_THREADS=16
srun ./ahmdal.exe

export OMP_NUM_THREADS=32
srun ./ahmdal.exe
