#!/usr/bin/env bash
#PBS -N Tife_gipaw
#PBS -o output_TiFe_gipaw.file
#PBS -e error_TiFe_gipaw.file
#PBS -l walltime=2:30:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=12gb

# Go to the directory where the job was submitted
cd "$PBS_O_WORKDIR"

# Load required module
module load QuantumESPRESSO


# Run gipaw
mpirun -np "$PBS_NP" gipaw.x -in TiFe.gipaw.in > TiFe.gipaw.out

