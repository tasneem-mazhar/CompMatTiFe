#!/bin/bash

# ============================================================
# Quantum ESPRESSO SCF automation with varying ecutwfc
# Author: Tasneem
# ============================================================

# Path to QE executable
PWX="/home/rafi/QE/qe-7.4.1/bin/pw.x"
# Pseudopotential directory
PSEUDO_DIR="/home/rafi/QE/qe-7.4.1/pseudo"
# MPI processes
NP=4
# Fixed K-points
K=10

# Create output directory if not exists
mkdir -p out

# Loop for 10 runs
for i in {1..15}; do
    PREFIX="TiFe${i}"
    INPUT="${PREFIX}.in"
    OUTPUT="${PREFIX}.out"

    # Define ecutwfc (e.g., 40, 50, 60, …)
    ECUTWFC=$((40 + (i-1)*10))   # 40,45,50,... if you want step=5
    ECUTRHO=$((ECUTWFC * 10))

    # Generate input file
    cat > "$INPUT" << EOF
&CONTROL
  calculation='scf',
  outdir='./out/',
  prefix='${PREFIX}',
  pseudo_dir='${PSEUDO_DIR}',
  verbosity='low',
  tprnfor=.true.,
  tstress=.true.,
/
&SYSTEM
  ibrav = 0,
  A = 2.93868,
  nat = 2,
  ntyp = 2,
  ecutwfc = ${ECUTWFC},
  ecutrho = ${ECUTRHO},
  input_dft = 'pbe',
  occupations = 'smearing',
  smearing = 'mv',
  degauss = 0.005d0,
/
&ELECTRONS
  conv_thr = 1d-08,
  mixing_beta = 0.7d0,
/
CELL_PARAMETERS {alat}
  1.000000000000000   0.000000000000000   0.000000000000000
  0.000000000000000   1.000000000000000   0.000000000000000
  0.000000000000000   0.000000000000000   1.000000000000000
ATOMIC_SPECIES
  Fe   55.84500  Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
  Ti   47.86700  ti_pbe_v1.4.uspp.F.UPF
ATOMIC_POSITIONS {crystal}
Fe   0.500000000000000   0.500000000000000   0.500000000000000
Ti   0.000000000000000   0.000000000000000   0.000000000000000
K_POINTS {automatic}
${K} ${K} ${K} 0 0 0
EOF

    # Run SCF calculation
    echo ">>> Running SCF for ${PREFIX} with ecutwfc=${ECUTWFC}, ecutrho=${ECUTRHO}, K=${K}x${K}x${K} ..."
    mpirun -np ${NP} ${PWX} -inp ${INPUT} > ${OUTPUT}
    echo ">>> Completed ${PREFIX}"
    echo "--------------------------------------------"
done

echo "All 10 ecutwfc convergence calculations finished!"
