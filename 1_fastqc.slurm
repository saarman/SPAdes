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

# Input folder(s) - must run for twice, once for each directory
input_folder=("/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702")         # also "UT-M70330-240718"
output_folder=("/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/fastqc_files")

# Create the output folder if it doesn't exist
mkdir -p ${output_folder}

# run fastqc for all .fastq.gz files in directory
fastqc *.fastq.gz -o ./fastqc_files

#change permissions
chmod -R g+w "$output_folder"



