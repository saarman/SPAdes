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

