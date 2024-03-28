#!/bin/bash
#SBATCH --account=class-cse4163
#SBATCH --qos=class-cse4163
#SBATCH --job-name=FluidOMP64
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --no-reque
OMP_NUM_THREADS=64 ./fluidomp -n 64 -o fkte64.dat >& out.64thread



