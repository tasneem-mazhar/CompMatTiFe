#!/bin/bash
# AUTOMATIC EOS GENERATION SCRIPT FOR QE

PW="/home/rafi/QE/qe-7.4.1/bin/pw.x"
NP=6

# Base lattice constant
A0=2.93868

# Percent changes for scaling (2%, 4%, … 20%)
percentages=(0.005 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.1)

# Output file for energy-volume data
EV_FILE="ev_data.dat"
echo "# Volume(Ang^3)    Energy(Ry)" > $EV_FILE

count=1

# =====================================================
# FUNCTION TO GENERATE INPUT FILES
# =====================================================
generate_input() {
Anew=$1
prefix=$2
infile=$3

cat > "$infile" << EOF
&CONTROL
  calculation='scf',
  outdir='./out',
  prefix='$prefix',
  pseudo_dir='/home/rafi/QE/qe-7.4.1/pseudo',
  verbosity='low',
  tprnfor=.true.,
  tstress=.true.,
/
&SYSTEM
  ibrav=0,
  A=$Anew,
  nat=2,
  ntyp=2,
  ecutwfc=140,
  ecutrho=1400,
  input_dft='pbe',
  occupations='smearing',
  smearing='mv',
  degauss=0.005d0,
/
&ELECTRONS
  conv_thr=1d-08,
  mixing_beta=0.7d0,
/
CELL_PARAMETERS {alat}
  1.0 0.0 0.0
  0.0 1.0 0.0
  0.0 0.0 1.0
ATOMIC_SPECIES
  Fe  55.845  Fe.pbe-spn-kjpaw_psl.0.2.1.UPF
  Ti  47.867  ti_pbe_v1.4.uspp.F.UPF
ATOMIC_POSITIONS {crystal}
Fe  0.5  0.5  0.5
Ti  0.0  0.0  0.0
K_POINTS {automatic}
10 10 10 0 0 0
EOF
}

# =====================================================
# RUN 10 EXPANDED (A increased)
# =====================================================
for p in "${percentages[@]}"; do
    Anew=$(echo "$A0*(1+$p)" | bc -l)
    infile="TiFe${count}.in"
    outfile="TiFe${count}.out"
    prefix="TiFe${count}"

    echo ">>> Generating $infile (A = $Anew)"
    generate_input "$Anew" "$prefix" "$infile"

    echo ">>> Running $prefix ..."
    mpirun -np $NP $PW -inp "$infile" > "$outfile"

    # Extract volume (Ang^3)
    VOL=$(grep -m 1 "unit-cell volume" "$outfile" | awk '{print $4}')
    # Extract energy (Ry)
    ENE=$(grep -m 1 "!    total energy" "$outfile" | awk '{print $5}')
# Extract Pressure (kbar)
    PRE=$(grep -m 1 "P=" "$outfile" | awk '{print $6}')
    echo "$VOL   $ENE   $PRE" >> "$EV_FILE"

    count=$((count+1))
done

# =====================================================
# RUN 10 COMPRESSED (A decreased)
# =====================================================
for p in "${percentages[@]}"; do
    Anew=$(echo "$A0*(1-$p)" | bc -l)
    infile="TiFe${count}.in"
    outfile="TiFe${count}.out"
    prefix="TiFe${count}"

    echo ">>> Generating $infile (A = $Anew)"
    generate_input "$Anew" "$prefix" "$infile"

    echo ">>> Running $prefix ..."
    mpirun -np $NP $PW -inp "$infile" > "$outfile"

    # Extract volume and energy
    VOL=$(grep -m 1 "unit-cell volume" "$outfile" | awk '{print $4}')
    ENE=$(grep -m 1 "!    total energy" "$outfile" | awk '{print $5}')

    echo "$VOL   $ENE" >> "$EV_FILE"

    count=$((count+1))
done

echo "=================================================="
echo " All calculations complete."
echo " Energy-volume data written to: $EV_FILE"
echo "=================================================="
