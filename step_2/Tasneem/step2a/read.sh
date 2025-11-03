#!/bin/bash

output_file="TiFe_vol_energy.txt"

# Overwrite/initialize output file (no header)
: > "$output_file"

for i in {1..20}; do
    file="TiFe${i}.out"

    if [[ -f "$file" ]]; then
        # Extract numeric unit-cell volume (a.u.^3)
        volume=$(grep -m1 "unit-cell volume" "$file" | sed -n 's/.*=\s*\([-0-9.eE+]\+\).*/\1/p')

        # Extract numeric total energy (Ry)
        energy=$(grep -m1 "! *total energy" "$file" | sed -n 's/.*=\s*\([-0-9.eE+]\+\).*/\1/p')

        # Fallback to NA if not found
        [[ -z "$volume" ]] && volume="NA"
        [[ -z "$energy" ]] && energy="NA"

        # Append volume and energy as two columns (tab-separated)
        printf "%s\t%s\n" "$volume" "$energy" >> "$output_file"
    else
        # File missing: warn to stderr (does not affect output file contents)
        echo "Warning: $file not found" >&2
    fi
done
