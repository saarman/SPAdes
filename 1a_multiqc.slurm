#!/bin/bash

#SBATCH --partition=saarman-np   
#SBATCH --account=saarman-np
#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=24576 # memory given in MB
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="multiqc"
#SBATCH --mail-user=emily.calhoun@usu.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Aggregate FastQC (v0.12.1) reports, run on 10k-read seqtk subsets (*_sub1)

# Load Module
module load multiqc/1.12

# Input and output files
in_path="/uufs/chpc.utah.edu/common/home/u1055819/saarman-group/uphlfiles/fastqc"
outdir="/uufs/chpc.utah.edu/common/home/u1055819/saarman-group/uphlfiles/fastqc/multiqc/"

# Create the output folder if it doesn't exist
mkdir -p "$outdir"

# Record version
multiqc --version > "${outdir}multiqc_version.txt"

multiqc "$in_path" --outdir "$outdir"

# Give the group write access to the output
chmod -R g+w "$outdir"
