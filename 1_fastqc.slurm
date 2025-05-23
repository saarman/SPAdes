#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --partition=saarman-shared-np   
#SBATCH --account=saarman-np
#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=24576 # memory given in MB
#SBATCH --nodes=1   # number of nodes
# #SBATCH --ntasks-per-node=16   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="fastqc"
# #SBATCH --mail-user=emily.calhoun@usu.edu   # insert email address
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

# Load Module
module load fastqc/0.12.1

# CD to working directory
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles

# Input folders
input_folders=("UT-M70330-240718" "UT-M07101-240702")

# Loop over each input folder
for input_folder in "${input_folders[@]}"; do
    echo "Processing folder: $input_folder"

    output_folder="${input_folder}/fastqfiles"
    mkdir -p "$output_folder"
    chmod -R g+w "$output_folder"

    # Use 'find' to get all .fastq or .fastq.gz files recursively
    find "$input_folder" -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | while read file; do
        echo "Running FastQC on $file"
        fastqc "$file" -o "$output_folder"
    done
done

#change permissions
chmod -R g+w "$output_folder"



