#!/bin/bash

# ============================================================
# Script to extract total energy and pressure from QE outputs
# Author: Tasneem
# ============================================================

OUTFILE="E_stress_summary.dat"

# Write header
echo -e "Run\tKpoints\tTotal_Energy(Ry)\tPressure(kbar)" > "$OUTFILE"

# Loop over output files
for i in {1..10}; do
    FILE="TiFe${i}.out"
    K=$((i * 2))

    if [[ -f "$FILE" ]]; then
        # Extract total energy (first occurrence)
        ENERGY=$(grep -m1 "!    total energy" "$FILE" | awk '{print $5}')
        
        # Extract pressure (P= value in kbar)
        PRESSURE=$(grep -m1 "total   stress" "$FILE" | awk -F'P=' '{print $2}' | awk '{print $1}')

        # Save to summary
        echo -e "${i}\t${K}x${K}x${K}\t${ENERGY}\t${PRESSURE}" >> "$OUTFILE"
    else
        echo "Warning: $FILE not found!"
    fi
done

echo "✅ Extraction complete! Data saved in ${OUTFILE}"
