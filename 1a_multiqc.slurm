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

base="/uufs/chpc.utah.edu/common/home/u1055819/saarman-group/uphlfiles"

for run in UT-M07101-240702 UT-M70330-240718; do
  in_path="${base}/${run}/fastqc"
  outdir="${base}/${run}/multiqc/"
  mkdir -p "$outdir"
  multiqc --version > "${outdir}multiqc_version.txt"
  multiqc "$in_path" --outdir "$outdir"
  chmod -R g+w "$outdir"
done
