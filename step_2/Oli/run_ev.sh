#!/bin/bash
# Simple E(V) curve script for TiFe
# Usage: bash run_ev.sh

# Check if TiFe.in exists
if [ ! -f "TiFe.in" ]; then
    echo "ERROR: TiFe.in not found in current directory!"
    exit 1
fi

# Configuration
SCALE_FACTORS="0.8 0.9 0.95 0.97 0.99 1.00 1.01 1.03 1.05 1.1 1.2"

# Get original lattice parameter from TiFe.in
ORIGINAL_A=$(grep "A =" TiFe.in | awk '{print $3}')
echo "Original lattice parameter: A = $ORIGINAL_A"

# Create directory and go there
mkdir -p ev_curve
cd ev_curve

# Create input files and job scripts
counter=1
for scale in $SCALE_FACTORS; do
    # Calculate new lattice parameter
    new_a=$(echo "$ORIGINAL_A * $scale" | bc -l)
    
    # Create input file by modifying original TiFe.in
    sed "s/A =.*/A = ${new_a}/" ../TiFe.in > input_${counter}.in
    sed -i "s/prefix='TiFe'/prefix='TiFe_${counter}'/" input_${counter}.in
    sed -i "s|pseudo_dir='\./|pseudo_dir='../|" input_${counter}.in
    
    # Create job script
    cat > job_${counter}.sh << EOF
#!/bin/sh
#PBS -N TiFe_${counter}
#PBS -o output_${counter}.file
#PBS -e error_${counter}.file
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=16
#PBS -l mem=4gb

cd \$PBS_O_WORKDIR
module load QuantumESPRESSO/7.4-foss-2024a
mpirun -np \$PBS_NP pw.x -input input_${counter}.in > output_${counter}.out
EOF

    chmod +x job_${counter}.sh
    
    # Submit job
    qsub job_${counter}.sh
    echo "Submitted job ${counter} with a=${new_a}"
    
    counter=$((counter + 1))
    sleep 0.5
done

echo ""
echo "All jobs submitted!"

