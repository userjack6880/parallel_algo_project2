#!/bin/bash
#SBATCH --account=class-cse4163
#SBATCH --qos=class-cse4163
#SBATCH --job-name=Fluid32t2p
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=32
#SBATCH --ntasks-per-node=2
#SBATCH --time=00:30:00
#SBATCH --no-reque
module load openmpi
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun -np 2 ./fluidmpi -n 64 -o fkte32t2p.dat >& out.mpi32t2p

