#!/bin/bash
#SBATCH --account=class-cse4163
#SBATCH --qos=class-cse4163
#SBATCH --job-name=Fluid2t32p
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:30:00
#SBATCH --no-reque
module load openmpi
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun -np 32 ./fluidmpi -n 64 -o fkte2t32p.dat >& out.mpi2t32p

