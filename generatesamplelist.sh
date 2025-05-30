#!/bin/bash

# Directories to scan for FASTQ files — update this list as needed
INPUT_DIRS=(
    "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702"
)

# Where to save the final list of sample names
SAMPLES_FILE="/uufs/chpc.utah.edu/common/home/saarman-group1/projects/bulkfastq/samples.txt"

# Make sure the output directory exists
mkdir -p "$(dirname "$SAMPLES_FILE")"

# Clear the file if it already exists
> "$SAMPLES_FILE"

# Loop through each input directory
for dir in "${INPUT_DIRS[@]}"; do
    if [[ -d "$dir" ]]; then
        echo "Scanning $dir ..."
        find "$dir" -type f -name "*_R1_001.fastq.gz" | while read -r f; do
            filename=$(basename "$f")
            # Get everything before the first underscore
            sample=$(echo "$filename" | cut -d_ -f1)
            echo "$sample" >> "$SAMPLES_FILE"
        done
    else
        echo "❌ Directory does not exist: $dir"
    fi
done

echo "✅ samples.txt written to: $SAMPLES_FILE"
wc -l "$SAMPLES_FILE"
