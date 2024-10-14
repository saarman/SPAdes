#!/bin/sh
#SBATCH --time=336:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20          # same as $max set in ForkManager
#SBATCH --account=saarman-np
#SBATCH --partition=saarman-shared-np   
#SBATCH --job-name=MMseqs2_try4
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=norah.saarman@usu.edu

# Load modules
module load mmseqs2/oct24  # change to module name

# Change to the directory where the input data is located
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input

# Concatenate all of the input into one file
rm all.fasta; cat *.fasta > all.fasta
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/

# Run the Perl script with the input files
perl /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/scripts/4a_MMseqs2.pl all.fasta

# Permissions
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/
