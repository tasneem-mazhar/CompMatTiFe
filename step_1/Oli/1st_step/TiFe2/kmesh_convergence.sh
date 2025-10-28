#!/usr/bin/env bash
#
# Job script for k-mesh convergence testing on HPC cluster
# Submit with: qsub kmesh_convergence_hpc.sh
#
#PBS -N TiFe_kmesh
#PBS -o kmesh_output.txt
#PBS -e kmesh_error.txt
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=8gb

# Change to directory where job was submitted
cd $PBS_O_WORKDIR

# Load Quantum Espresso module
module load QuantumESPRESSO/7.4-foss-2024a

# Print job information
echo "Job started on: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $PBS_JOBID"
echo "Working directory: $PBS_O_WORKDIR"
echo ""

# Base input file
BASE_INPUT="TiFe.in"

# Create results directory
RESULTS_DIR="TiFe_kmesh_convergence"
mkdir -p "$RESULTS_DIR"

echo "════════════════════════════════════════════════════════════"
echo "  K-mesh Convergence Test (HPC)"
echo "════════════════════════════════════════════════════════════"
echo ""

# K-mesh values to test
K_MESHES=(
    "2 2 2"
    "3 3 3"
    "4 4 4"
    "5 5 5"
    "6 6 6"
    "7 7 7"
    "8 8 8"
    "9 9 9"
)

# Results file
RESULTS_FILE="${RESULTS_DIR}/kmesh_convergence.dat"
echo "# k-mesh    Pressure(kbar)    Energy(Ry)    Time(s)" > "$RESULTS_FILE"

# Process each k-mesh
for kmesh in "${K_MESHES[@]}"; do
    read k1 k2 k3 <<< "$kmesh"
    mesh_label="${k1}x${k2}x${k3}"
    
    echo "Testing k-mesh: $mesh_label"
    
    INPUT_FILE="${RESULTS_DIR}/TiFe2_k${k1}${k2}${k3}.in"
    OUTPUT_FILE="${INPUT_FILE%.in}.out"
    
    # Modify k-mesh in input file
    awk -v k1="$k1" -v k2="$k2" -v k3="$k3" '
    BEGIN { in_kpoints = 0; kpoints_written = 0 }
    {
        if ($0 ~ /^K_POINTS/) {
            in_kpoints = 1
            print $0
            next
        }
        
        if (in_kpoints && !kpoints_written) {
            print k1 " " k2 " " k3 " 0 0 0"
            kpoints_written = 1
            in_kpoints = 0
            next
        }
        
        if (in_kpoints && kpoints_written) {
            in_kpoints = 0
            next
        }
        
        print $0
    }
    ' "$BASE_INPUT" > "$INPUT_FILE"
    
    # Run calculation with MPI
    START_TIME=$(date +%s)
    
    mpirun -np $PBS_NP pw.x -input "$INPUT_FILE" > "$OUTPUT_FILE" 2>&1
    
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    
    # Extract results
    PRESSURE=$(grep -a "total   stress" "$OUTPUT_FILE" -A 3 | grep "P=" | awk '{print $6}')
    ENERGY=$(grep -a "!    total energy" "$OUTPUT_FILE" | tail -1 | awk '{print $5}')
    
    if [ -z "$PRESSURE" ]; then
        PRESSURE="N/A"
    fi
    
    if [ -z "$ENERGY" ]; then
        ENERGY="N/A"
    fi
    
    echo "  ${mesh_label}: P = ${PRESSURE} kbar, E = ${ENERGY} Ry, Time = ${ELAPSED}s"
    echo "${mesh_label}    ${PRESSURE}    ${ENERGY}    ${ELAPSED}" >> "$RESULTS_FILE"
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
