#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --partition=saarman-shared-np   
#SBATCH --account=saarman-np
#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=24576 # memory given in MB
#SBATCH --nodes=1   # number of nodes
# #SBATCH --ntasks-per-node=16   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="SPAdes"
# #SBATCH --mail-user=emily.calhoun@usu.edu   # email address
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

# Load SPAdes
module load spades

# Define paths
SAMPLES_FILE="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdesscripts/samples.txt"
OUTPUT_BASE="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdesscripts/denovo_assembly"

# Get this task's line from samples.txt
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES_FILE")
INPUT_DIR=$(echo "$LINE" | awk '{print $1}')
SAMPLE=$(echo "$LINE" | awk '{print $2}')
OUTPUT_DIR="$OUTPUT_BASE/${SAMPLE}"
mkdir -p "$OUTPUT_DIR"

# Skip if output already exists
if [[ -f "$OUTPUT_DIR/contigs.fasta" ]]; then
  echo "[$SAMPLE] Already complete. Skipping."
  exit 0
fi

# Run SPAdes
echo "[$SAMPLE] Running SPAdes..."
cd "$INPUT_DIR" || exit 1

spades.py \
  -1 "${SAMPLE}_L001_R1_001.fastq.gz" \
  -2 "${SAMPLE}_L001_R2_001.fastq.gz" \
  -o "$OUTPUT_DIR" \
  --isolate \
  -t 8 -m 32

chmod -R g+w "$OUTPUT_DIR"
echo "[$SAMPLE] Done."



