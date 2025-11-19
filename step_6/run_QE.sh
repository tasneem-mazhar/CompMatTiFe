#!/usr/bin/env bash
#PBS -N Tife_magnetism_vcrelax
#PBS -o output_TiFe_HH.file
#PBS -e error_TiFe_HH.file
#PBS -l walltime=2:30:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=12gb

# Go to the directory where the job was submitted
cd "$PBS_O_WORKDIR"

# Load required module
module load QuantumESPRESSO


# Run Quantum ESPRESSO for each input file
mpirun -np "$PBS_NP" pw.x -input TiFe_optimised.in > TiFe_optimised_magn.out
