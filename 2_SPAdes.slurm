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

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load spades

# Define array of input folders
input_folders=(
  "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M70330-240718"
  "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702"
)

output_folder="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/denovo_assembly"

mkdir -p "$output_folder"

for input_folder in "${input_folders[@]}"; do
  cd "$input_folder" || { echo "Cannot enter $input_folder"; continue; }
  
  for SAMPLE in $(ls -1 *R1_001.fastq.gz | cut -d_ -f1-2); do
    echo "Running SPAdes for sample: $SAMPLE in $input_folder"
    spades.py \
      -1 "${SAMPLE}_L001_R1_001.fastq.gz" \
      -2 "${SAMPLE}_L001_R2_001.fastq.gz" \
      -o "${output_folder}/${SAMPLE}" \
      --isolate
  done
done

chmod -R g+w "$output_folder"




