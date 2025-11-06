#!/bin/bash
# Extract E(V) data from QE output files
# Run this inside the ev_curve directory after jobs finish

cd ev_curve
echo "# Volume(Ang^3)  Energy(Ry)" > ev_data.txt

for i in {1..11}; do
    outfile="output_${i}.out"
    
    if [ -f "$outfile" ] && grep -q "JOB DONE" "$outfile"; then
        volume=$(grep "unit-cell volume" "$outfile" | tail -1 | awk '{print $4}')
        energy=$(grep "!" "$outfile" | tail -1 | awk '{print $5}')
        
        if [ -n "$volume" ] && [ -n "$energy" ]; then
            printf "%14.6f  %18.10f\n" $volume $energy >> ev_data.txt
            echo "Point $i: V=$volume Ang^3, E=$energy Ry"
        fi
    else
        echo "Point $i: not finished or failed"
    fi
done

echo ""
echo "Data saved to: ev_curve/ev_data.txt"
cat ev_data.txt
