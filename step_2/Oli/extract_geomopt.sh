#!/bin/bash
# Extract results from two-step geometry optimization
# Usage: bash extract_geomopt.sh

echo "======================================================================"
echo "GEOMETRY OPTIMIZATION RESULTS"
echo "======================================================================"
echo ""

# ======================================================================
# STEP 1: RELAX RESULTS
# ======================================================================

echo "STEP 1: RELAX (Position Optimization)"
echo "----------------------------------------------------------------------"

if [ -f "TiFe_relax.out" ]; then
    if grep -q "JOB DONE" "TiFe_relax.out"; then
        echo "Status: COMPLETED"
        echo ""
        grep "!" "TiFe_relax.out" | tail -1
        echo ""
        echo "Final forces:"
        awk '/Begin final coordinates/,/End final coordinates/' "TiFe_relax.out" | grep -A 5 "Forces acting" | head -6
    else
        echo "Status: NOT COMPLETED or FAILED"
    fi
else
    echo "File not found: TiFe_relax.out"
fi

echo ""
echo "======================================================================"
echo ""

# ======================================================================
# STEP 2: VC-RELAX RESULTS
# ======================================================================

echo "STEP 2: VC-RELAX (Full Optimization)"
echo "----------------------------------------------------------------------"

outfile="TiFe_vcrelax.out"

if [ ! -f "$outfile" ]; then
    echo "File not found: $outfile"
    echo "Step 2 may not have run yet or is still running."
    echo ""
    exit 1
fi

if ! grep -q "JOB DONE" "$outfile"; then
    echo "WARNING: Calculation may not have finished successfully!"
    echo ""
fi

echo "Final Energy:"
echo "----------------------------------------------------------------------"
grep "!" "$outfile" | tail -1
echo ""

echo "Final Volume:"
echo "----------------------------------------------------------------------"
grep "unit-cell volume" "$outfile" | tail -1
volume=$(grep "unit-cell volume" "$outfile" | tail -1 | awk '{print $4}')
echo ""

echo "Final Lattice Parameters:"
echo "----------------------------------------------------------------------"
alat=$(grep "lattice parameter (alat)" "$outfile" | tail -1 | awk '{print $5}')
echo "alat = $alat bohr"
a_ang=$(echo "$alat * 0.529177" | bc -l)
echo "a = $a_ang Angstrom (for cubic cell)"
echo ""

echo "CELL_PARAMETERS:"
awk '/Begin final coordinates/,/End final coordinates/' "$outfile" | grep -A 3 "CELL_PARAMETERS"
echo ""

echo "ATOMIC_POSITIONS:"
awk '/Begin final coordinates/,/End final coordinates/' "$outfile" | grep -A 10 "ATOMIC_POSITIONS"
echo ""

echo "Final Forces:"
echo "----------------------------------------------------------------------"
awk '/Begin final coordinates/,/End final coordinates/' "$outfile" | grep -A 5 "Forces acting" | head -6
echo ""

echo "Final Stress Tensor:"
echo "----------------------------------------------------------------------"
grep -A 3 "total.*stress" "$outfile" | tail -4
echo ""

echo "======================================================================"
echo "SUMMARY"
echo "======================================================================"
echo "Optimized volume: $volume Ang^3"
echo "Optimized lattice parameter: $a_ang Angstrom"
echo ""
echo "Compare this volume with V0 from your E(V) fit!"
echo "======================================================================"

# Save to file
{
    echo "VC-RELAX Results for TiFe"
    echo "=========================="
    echo ""
    echo "Final volume: $volume Ang^3"
    echo "Final lattice parameter: $a_ang Angstrom"
    echo ""
    grep "!" "$outfile" | tail -1
} > geomopt_results.txt

echo ""
echo "Results saved to: geomopt_results.txt"
