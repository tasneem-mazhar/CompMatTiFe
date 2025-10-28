#!/usr/bin/env bash
#
# Job script for ecutrho factor convergence testing on HPC cluster
# IMPORTANT: Edit KMESH and ECUTWFC variables below with converged values
# Submit with: qsub ecutrho_convergence_hpc.sh
#
#PBS -N TiFe2_ecutrho
#PBS -o ecutrho_output.txt
#PBS -e ecutrho_error.txt
#PBS -l walltime=4:00:00
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

# EDIT THESE: Use converged values from previous tests
KMESH="9 9 9"
ECUTWFC=140

BASE_INPUT="TiFe.in"
RESULTS_DIR="TiFe_ecutrho_convergence"
mkdir -p "$RESULTS_DIR"

echo "════════════════════════════════════════════════════════════"
echo "  ecutrho Factor Convergence Test (HPC)"
echo "════════════════════════════════════════════════════════════"
echo ""
echo "Using converged k-mesh: $KMESH"
echo "Using converged ecutwfc: $ECUTWFC Ry"
echo ""

# Factors to test
FACTORS=(2 4 5 6 7 8 9 10 12 15 20 25 30)

# Results file
RESULTS_FILE="${RESULTS_DIR}/ecutrho_factor_convergence.dat"
echo "# Factor    ecutrho(Ry)    Pressure(kbar)    Energy(Ry)    Time(s)" > "$RESULTS_FILE"

# Process each factor
for factor in "${FACTORS[@]}"; do
    ecutrho=$((ECUTWFC * factor))
    
    echo "Testing factor = $factor (ecutrho = $ecutrho Ry)"
    
    INPUT_FILE="${RESULTS_DIR}/TiFe2_factor${factor}.in"
    OUTPUT_FILE="${INPUT_FILE%.in}.out"
    
    # Modify input file
    # Parse KMESH into components
    read k1 k2 k3 <<< "$KMESH"
    
    awk -v ecutwfc="$ECUTWFC" -v ecutrho="$ecutrho" -v k1="$k1" -v k2="$k2" -v k3="$k3" '
    BEGIN { 
        in_system = 0
        in_kpoints = 0
        kpoints_done = 0
    }
    {
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
        
        if (in_system && $0 ~ /ecutwfc/) {
            print "   ecutwfc=" ecutwfc ","
            next
        }
        
        if (in_system && $0 ~ /ecutrho/) {
            print "   ecutrho=" ecutrho ","
            next
        }
        
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
        
        print $0
    }
    ' "$BASE_INPUT" > "$INPUT_FILE"
    
    # Verify file was created
    if [ ! -s "$INPUT_FILE" ]; then
        echo "  ERROR: Failed to create input file!"
        continue
    fi
    
    # Run calculation
    START_TIME=$(date +%s)
    
    echo "  Running: mpirun -np $PBS_NP pw.x -input $INPUT_FILE"
    mpirun -np $PBS_NP pw.x -input "$INPUT_FILE" > "$OUTPUT_FILE" 2>&1
    EXIT_CODE=$?
    
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    
    echo "  Exit code: $EXIT_CODE, Time: ${ELAPSED}s"
    
    # Check output file
    if [ ! -s "$OUTPUT_FILE" ]; then
        echo "  ERROR: Output file is empty!"
        PRESSURE="EMPTY"
        ENERGY="EMPTY"
    else
        # Extract results
        PRESSURE=$(grep -a "total   stress" "$OUTPUT_FILE" -A 3 | grep "P=" | awk '{print $6}')
        ENERGY=$(grep -a "!    total energy" "$OUTPUT_FILE" | tail -1 | awk '{print $5}')
        
        if [ -z "$PRESSURE" ]; then
            PRESSURE="N/A"
        fi
        
        if [ -z "$ENERGY" ]; then
            ENERGY="N/A"
        fi
    fi
    
    echo "  P = ${PRESSURE} kbar, E = ${ENERGY} Ry, Time = ${ELAPSED}s"
    echo "${factor}    ${ecutrho}    ${PRESSURE}    ${ENERGY}    ${ELAPSED}" >> "$RESULTS_FILE"
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
