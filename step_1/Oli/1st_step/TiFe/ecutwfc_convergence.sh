#!/usr/bin/env bash
#
# Job script for ecutwfc convergence testing on HPC cluster
# IMPORTANT: Edit KMESH variable below with converged value from k-mesh test
# Submit with: qsub ecutwfc_convergence_hpc.sh
#
#PBS -N TiFe2_ecutwfc
#PBS -o ecutwfc_output.txt
#PBS -e ecutwfc_error.txt
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=16
#PBS -l mem=16gb

# Change to directory where job was submitted
cd $PBS_O_WORKDIR

# Load Quantum Espresso module
module load QuantumESPRESSO/7.4-foss-2024a

# Print job information
echo "Job started on: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $PBS_JOBID"
echo ""

# EDIT THIS: Use converged k-mesh from previous test
KMESH="9 9 9"

BASE_INPUT="TiFe.in"
RESULTS_DIR="TiFe_ecutwfc_convergence"
mkdir -p "$RESULTS_DIR"

echo "════════════════════════════════════════════════════════════"
echo "  ecutwfc Convergence Test (HPC)"
echo "════════════════════════════════════════════════════════════"
echo ""
echo "Using converged k-mesh: $KMESH"
echo ""

# ecutwfc values to test (ecutrho = 5 × ecutwfc)
ECUTWFC_VALUES=(20 40 50 60 70 80 100 120 140 170 200 250)

# Results file
RESULTS_FILE="${RESULTS_DIR}/ecutwfc_convergence.dat"
echo "# ecutwfc(Ry)    ecutrho(Ry)    Pressure(kbar)    Energy(Ry)    Time(s)" > "$RESULTS_FILE"

# Process each ecutwfc value
for ecutwfc in "${ECUTWFC_VALUES[@]}"; do
    ecutrho=$((ecutwfc * 5))
    
    echo "Testing ecutwfc = $ecutwfc Ry (ecutrho = $ecutrho Ry)"
    
    INPUT_FILE="${RESULTS_DIR}/TiFe2_ecut${ecutwfc}.in"
    OUTPUT_FILE="${INPUT_FILE%.in}.out"
    
    # Modify input file
    # Parse KMESH into components
    read k1 k2 k3 <<< "$KMESH"
    
    awk -v ecutwfc="$ecutwfc" -v ecutrho="$ecutrho" -v k1="$k1" -v k2="$k2" -v k3="$k3" '
    BEGIN { 
        in_system = 0
        in_kpoints = 0
        kpoints_done = 0
    }
    {
        # Handle SYSTEM block
        if ($0 ~ /^&SYSTEM/) {
            in_system = 1
            print $0
            next
        }
        
        if (in_system && $0 ~ /^\//) {
            in_system = 0
            print $0
            next
        }
        
        # Replace ecutwfc in SYSTEM block
        if (in_system && $0 ~ /ecutwfc/) {
            print "   ecutwfc=" ecutwfc ","
            next
        }
        
        # Replace ecutrho in SYSTEM block
        if (in_system && $0 ~ /ecutrho/) {
            print "   ecutrho=" ecutrho ","
            next
        }
        
        # Handle K_POINTS
        if ($0 ~ /^K_POINTS/) {
            in_kpoints = 1
            print $0
            next
        }
        
        if (in_kpoints && !kpoints_done) {
            print k1 " " k2 " " k3 " 0 0 0"
            kpoints_done = 1
            in_kpoints = 0
            next
        }
        
        # Skip original k-points line if already written
        if (in_kpoints && kpoints_done) {
            in_kpoints = 0
            next
        }
        
        # Default: print line as-is
        print $0
    }
    ' "$BASE_INPUT" > "$INPUT_FILE"
    
    # Verify input file was created and has content
    if [ ! -s "$INPUT_FILE" ]; then
        echo "  ERROR: Failed to create input file or file is empty!"
        continue
    fi
    
    # Run calculation
    START_TIME=$(date +%s)
    
    echo "  Running: mpirun -np $PBS_NP pw.x -input $INPUT_FILE > $OUTPUT_FILE"
    mpirun -np $PBS_NP pw.x -input "$INPUT_FILE" > "$OUTPUT_FILE" 2>&1
    EXIT_CODE=$?
    
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    
    echo "  Calculation finished with exit code: $EXIT_CODE, Time: ${ELAPSED}s"
    
    # Check if output file exists and has content
    if [ ! -s "$OUTPUT_FILE" ]; then
        echo "  ERROR: Output file is empty or does not exist!"
        echo "${ecutwfc}    ${ecutrho}    EMPTY_OUTPUT    EMPTY_OUTPUT    ${ELAPSED}" >> "$RESULTS_FILE"
        continue
    fi
    
    echo "  Output file size: $(wc -l < "$OUTPUT_FILE") lines"
    
    # Extract results
    PRESSURE=$(grep -a "total   stress" "$OUTPUT_FILE" -A 3 | grep "P=" | awk '{print $6}')
    ENERGY=$(grep -a "!    total energy" "$OUTPUT_FILE" | tail -1 | awk '{print $5}')
    
    if [ -z "$PRESSURE" ]; then
        PRESSURE="N/A"
    fi
    
    if [ -z "$ENERGY" ]; then
        ENERGY="N/A"
    fi
    
    echo "  P = ${PRESSURE} kbar, E = ${ENERGY} Ry, Time = ${ELAPSED}s"
    echo "${ecutwfc}    ${ecutrho}    ${PRESSURE}    ${ENERGY}    ${ELAPSED}" >> "$RESULTS_FILE"
done

echo ""
echo "════════════════════════════════════════════════════════════"
echo "  RESULTS SUMMARY"
echo "════════════════════════════════════════════════════════════"
echo ""
cat "$RESULTS_FILE"

echo ""
echo "Job completed on: $(date)"
echo "Results saved to: $RESULTS_FILE"
