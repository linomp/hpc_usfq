#!/bin/bash
#
#SBATCH --job-name=cuda_[tu_nombre]
#SBATCH --output=res.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 16
#SBATCH --nodelist compute-0-12
#SBATCH --time=1:00


srun gcc -std=c99 -fopenmp evaluar_func.c -o ahmdal.exe -lm

export OMP_NUM_THREADS=1
srun ./ahmdal.exe
