#!/bin/bash

# Output file with sample list
SAMPLES_FILE="samples.txt"
> "$SAMPLES_FILE"

# Input folders with FASTQ files
INPUT_DIRS=(
  "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M70330-240718"
  "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702"
)

for folder in "${INPUT_DIRS[@]}"; do
  echo "Scanning $folder"
  cd "$folder" || continue
  for r1 in *_R1_001.fastq.gz; do
    [[ -e "$r1" ]] || continue
    SAMPLE=$(echo "$r1" | sed 's/_L001_R1_001.fastq.gz//')
    echo "$folder $SAMPLE" >> "$SAMPLES_FILE"
  done
done

echo "✅ Sample list written to $SAMPLES_FILE"
