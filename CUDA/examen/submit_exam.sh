#!/bin/bash
#
#SBATCH --job-name=cuda_examen_lmp
#SBATCH --output=res.txt
#
#SBATCH --ntasks=1
#SBATCH --nodelist compute-0-12
#SBATCH --time=1:00

export PATH=$PATH:/usr/local/cuda/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib

srun nvcc -arch=sm_35 suma_vec.cu -o sol_cuda
srun ./sol_cuda
