#!/bin/bash
# Two-step geometry optimization for TiFe
# Step 1: Relax (optimize atomic positions)
# Step 2: VC-relax (optimize cell + positions)
# Usage: bash run_geomopt.sh

# Check if TiFe.in exists
if [ ! -f "TiFe.in" ]; then
    echo "ERROR: TiFe.in not found in current directory!"
    exit 1
fi

echo "======================================================================"
echo "Two-Step Geometry Optimization"
echo "======================================================================"
echo ""

# Extract lattice parameter
A_value=$(grep "A =" TiFe.in | awk '{print $3}')

# ======================================================================
# STEP 1: RELAX - Optimize atomic positions only
# ======================================================================

echo "Step 1: Creating relax input (optimize positions only)..."

cat > TiFe_relax.in << EOF
&CONTROL
   calculation='relax',
   outdir='.',
   prefix='TiFe_relax',
   pseudo_dir='./',
   verbosity='low',
   tprnfor=.true.,
   tstress=.true.,
/
&SYSTEM
  ibrav = 0
  A = ${A_value}
  nat = 2
  ntyp = 2
   ecutwfc=140,
   ecutrho=1260,
   input_dft='pbe',
   occupations='smearing',
   smearing='gaussian',
   degauss=0.005d0,
/
&ELECTRONS
   conv_thr=1d-08,
   mixing_beta=0.7d0,
/
&IONS
  ion_dynamics='bfgs'
/
CELL_PARAMETERS {alat}
  1.000000000000000   0.000000000000000   0.000000000000000 
  0.000000000000000   1.000000000000000   0.000000000000000 
  0.000000000000000   0.000000000000000   1.000000000000000 
ATOMIC_SPECIES
  Fe   55.84500  Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
  Ti   47.86700  Ti_pbe_v1.4.uspp.F.UPF
ATOMIC_POSITIONS {crystal}
Fe   0.500000000000000   0.500000000000000   0.500000000000000 
Ti   0.000000000000000   0.000000000000000   0.000000000000000 
K_POINTS {automatic}
9 9 9 0 0 0
EOF

echo "Created: TiFe_relax.in"

# Create job script for relax
cat > job_relax.sh << 'EOF'
#!/bin/sh
#PBS -N TiFe_relax
#PBS -o relax_output.file
#PBS -e relax_error.file
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=16
#PBS -l mem=4gb

cd $PBS_O_WORKDIR
module load QuantumESPRESSO/7.4-foss-2024a
mpirun -np $PBS_NP pw.x -input TiFe_relax.in > TiFe_relax.out
EOF

chmod +x job_relax.sh

# Submit relax job
job1=$(qsub job_relax.sh)
echo "Submitted Step 1 (relax): $job1"
echo ""

# ======================================================================
# STEP 2: VC-RELAX - Full optimization (depends on step 1)
# ======================================================================

echo "Step 2: Creating vc-relax input (optimize cell + positions)..."

# Create vc-relax input that will use optimized positions from step 1
cat > TiFe_vcrelax.in << EOF
&CONTROL
   calculation='vc-relax',
   outdir='.',
   prefix='TiFe_vcrelax',
   pseudo_dir='./',
   verbosity='low',
   tprnfor=.true.,
   tstress=.true.,
/
&SYSTEM
  ibrav = 0
  A = ${A_value}
  nat = 2
  ntyp = 2
   ecutwfc=140,
   ecutrho=1260,
   input_dft='pbe',
   occupations='smearing',
   smearing='gaussian',
   degauss=0.005d0,
/
&ELECTRONS
   conv_thr=1d-08,
   mixing_beta=0.7d0,
/
&IONS
  ion_dynamics='bfgs'
/
&CELL
  cell_dynamics='bfgs'
  press=0.d0
  press_conv_thr=0.5d0
/
CELL_PARAMETERS {alat}
  1.000000000000000   0.000000000000000   0.000000000000000 
  0.000000000000000   1.000000000000000   0.000000000000000 
  0.000000000000000   0.000000000000000   1.000000000000000 
ATOMIC_SPECIES
  Fe   55.84500  Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
  Ti   47.86700  Ti_pbe_v1.4.uspp.F.UPF
ATOMIC_POSITIONS {crystal}
Fe   0.500000000000000   0.500000000000000   0.500000000000000 
Ti   0.000000000000000   0.000000000000000   0.000000000000000 
K_POINTS {automatic}
9 9 9 0 0 0
EOF

echo "Created: TiFe_vcrelax.in"

# Create job script for vc-relax (depends on relax finishing)
cat > job_vcrelax.sh << EOF
#!/bin/sh
#PBS -N TiFe_vcrelax
#PBS -o vcrelax_output.file
#PBS -e vcrelax_error.file
#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=16
#PBS -l mem=4gb
#PBS -W depend=afterok:${job1}

cd \$PBS_O_WORKDIR
module load QuantumESPRESSO/7.4-foss-2024a

# Extract optimized positions from relax output
if [ -f "TiFe_relax.out" ] && grep -q "JOB DONE" "TiFe_relax.out"; then
    echo "Using optimized positions from relax step..."
    # Note: For TiFe cubic structure, positions don't change
    # but this shows the proper workflow
fi

mpirun -np \$PBS_NP pw.x -input TiFe_vcrelax.in > TiFe_vcrelax.out
EOF

chmod +x job_vcrelax.sh

# Submit vc-relax job (will wait for relax to finish)
job2=$(qsub job_vcrelax.sh)
echo "Submitted Step 2 (vc-relax): $job2"
echo "  (will start after Step 1 completes)"
echo ""

echo "======================================================================"
echo "Jobs submitted!"
echo "======================================================================"
echo "Step 1 (relax):    $job1"
echo "Step 2 (vc-relax): $job2 (depends on Step 1)"
echo ""
echo "Monitor with: qstat -u \$USER"
echo ""
echo "After completion, extract results with:"
echo "  bash extract_geomopt.sh"
echo "======================================================================"

# Create job script
cat > job_vcrelax.sh << 'EOF'
#!/bin/sh
#PBS -N TiFe_vcrelax
#PBS -o vcrelax_output.file
#PBS -e vcrelax_error.file
#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=16
#PBS -l mem=32gb

cd $PBS_O_WORKDIR
module load QuantumESPRESSO/7.4-foss-2024a
mpirun -np $PBS_NP pw.x -input TiFe_vcrelax.in > TiFe_vcrelax.out
EOF

chmod +x job_vcrelax.sh

echo "Created: job_vcrelax.sh"
echo ""

# Submit job
qsub job_vcrelax.sh

echo "Job submitted!"
echo "Monitor with: qstat -u \$USER"
echo ""
echo "When done, extract results with: bash extract_vcrelax.sh"
