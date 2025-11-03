#!/bin/bash

# Output summary file
output_file="info.txt"

# Header
echo -e "Run\tTotal_Energy(Ry)\tVolume(a.u.^3)\tPressure(kbar)" > "$output_file"

# Loop through all files TiFe1.out ... TiFe20.out
for i in {1..20}; do
    file="TiFe${i}.out"

    if [[ -f "$file" ]]; then
        # Extract total energy (Ry)
        energy=$(grep "! *total energy" "$file" | awk '{for(j=1;j<=NF;j++) if ($j=="=") print $(j+1)}')

        # Extract unit-cell volume (a.u.^3)
        volume=$(grep "unit-cell volume" "$file" | awk '{for(j=1;j<=NF;j++) if ($j=="=") print $(j+1)}')

        # Extract pressure (kbar)
        pressure=$(grep "(kbar)" "$file" | sed -n 's/.*P=\s*\([-0-9.]*\).*/\1/p')

        # Append to output
        echo -e "${i}\t${energy}\t${volume}\t${pressure}" >> "$output_file"
    else
        echo "Warning: $file not found" >&2
    fi
done

echo "✅ Extraction complete. Results saved in $output_file"
