#!/bin/bash

# Define the paths to your input FASTQ folders
INPUT_DIRS=(
  "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702"
  "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M06363-240702"
)

# Path where samples.txt will be saved
SAMPLES_FILE="/uufs/chpc.utah.edu/common/home/saarman-group1/your_project_dir/samples.txt"
> "$SAMPLES_FILE"  # Clear file

# Loop over input directories
for folder in "${INPUT_DIRS[@]}"; do
  cd "$folder" || continue

  for r1 in *_R1_001.fastq.gz; do
    [[ -e "$r1" ]] || continue
    sample=$(echo "$r1" | sed 's/_L001_R1_001.fastq.gz//')
    echo "$folder $sample" >> "$SAMPLES_FILE"
  done
done

echo "✅ samples.txt written to: $SAMPLES_FILE"
wc -l "$SAMPLES_FILE"

