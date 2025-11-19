#!/usr/bin/env bash
#PBS -N Tife_scf
#PBS -o output_TiFe_scf.file
#PBS -e error_TiFe_scf.file
#PBS -l walltime=2:30:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=12gb

# Go to the directory where the job was submitted
cd "$PBS_O_WORKDIR"

# Load required module
module load QuantumESPRESSO


# Run Quantum ESPRESSO for each input file
mpirun -np "$PBS_NP" pw.x -input TiFe.scf.in > TiFe.scf.out
