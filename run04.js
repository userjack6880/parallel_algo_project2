#!/bin/bash
#SBATCH --account=class-cse4163
#SBATCH --qos=class-cse4163
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --no-reque
OMP_NUM_THREADS=4 ./fluidomp -n 64 -o fkte04.dat >& data/run$1/out.04thread



