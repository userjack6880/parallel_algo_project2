#!/bin/bash
#SBATCH --account=class-cse4163
#SBATCH --qos=class-cse4163
#SBATCH --job-name=Fluid4t16p
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:30:00
#SBATCH --no-reque
module load openmpi
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun -np 16 ./fluidmpi -n 64 -o fkte4t16p.dat >& out.mpi4t16p

