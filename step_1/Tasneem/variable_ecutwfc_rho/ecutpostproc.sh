#!/bin/bash

# ============================================================
# Extract total energy and pressure from QE outputs
# for ecutwfc convergence study
# ============================================================

OUTFILE="E_pressure_ecut.dat"

# Header
echo -e "Run\tecutwfc\tecutrho\tTotal_Energy(Ry)\tPressure(kbar)" > "$OUTFILE"

# Loop over 10 runs
for i in {1..15}; do
    FILE="TiFe${i}.out"

    # ecutwfc and ecutrho must match what was used
    ECUTWFC=$((40 + (i-1)*10))
    ECUTRHO=$((ECUTWFC * 10))

    if [[ -f "$FILE" ]]; then
        # Extract total energy
        ENERGY=$(grep -m1 "!    total energy" "$FILE" | awk '{print $5}')
        # Extract pressure (P= value)
        PRESSURE=$(grep -m1 "total   stress" "$FILE" | awk -F'P=' '{print $2}' | awk '{print $1}')

        # Save to summary
        echo -e "${i}\t${ECUTWFC}\t${ECUTRHO}\t${ENERGY}\t${PRESSURE}" >> "$OUTFILE"
    else
        echo "Warning: $FILE not found!"
    fi
done

echo "✅ Extraction complete! Data saved in ${OUTFILE}"
