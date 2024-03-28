#!/bin/bash
#SBATCH --account=class-cse4163
#SBATCH --qos=class-cse4163
#SBATCH --job-name=Fluid8t8p
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:30:00
#SBATCH --no-reque
module load openmpi
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun -np 8 ./fluidmpi -n 64 -o fkte8t8p.dat >& out.mpi8t8p

