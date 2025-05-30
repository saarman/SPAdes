#!/bin/bash

# Set your project directory here (update this to your real project folder name)
PROJECT_DIR="/uufs/chpc.utah.edu/common/home/saarman-group1/YOUR_PROJECT_DIR"
INPUT_DIR_1="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M70330-240718"
INPUT_DIR_2="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702"

SAMPLES_FILE="$PROJECT_DIR/samples.txt"

# Make sure the project directory exists
mkdir -p "$PROJECT_DIR"

# Clear the samples.txt file
> "$SAMPLES_FILE"

# Function to process input directory
process_dir() {
  local INPUT_DIR="$1"
  if [ ! -d "$INPUT_DIR" ]; then
    echo "❌ ERROR: Directory does not exist: $INPUT_DIR"
    return
  fi

  echo "🔍 Processing directory: $INPUT_DIR"
  for R1 in "$INPUT_DIR"/*_R1_001.fastq.gz; do
    if [ -f "$R1" ]; then
      R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
      if [ -f "$R2" ]; then
        SAMPLE_NAME=$(basename "$R1" | cut -d'_' -f1)
        echo -e "${SAMPLE_NAME}\t${R1}\t${R2}" >> "$SAMPLES_FILE"
      else
        echo "⚠️ WARNING: R2 file missing for $R1"
      fi
    fi
  done
}

# Process both directories
process_dir "$INPUT_DIR_1"
process_dir "$INPUT_DIR_2"

echo "✅ samples.txt written to: $SAMPLES_FILE"
wc -l "$SAMPLES_FILE"
