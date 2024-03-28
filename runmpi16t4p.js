#!/bin/bash
#SBATCH --account=class-cse4163
#SBATCH --qos=class-cse4163
#SBATCH --job-name=Fluid16t4p
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --no-reque
module load openmpi
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun -np 4 ./fluidmpi -n 64 -o fkte16t4p.dat >& out.mpi16t4p

