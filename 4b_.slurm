#!/bin/sh
#SBATCH --time=24:00:00                # Set the time limit for the job
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks=1                     # Number of tasks
#SBATCH --cpus-per-task=4              # Number of CPU cores per task (adjust based on your needs)
#SBATCH --mem=128G                      # Memory allocation (adjust based on your needs)
#SBATCH --account=saarman-np           # Your account
#SBATCH --partition=saarman-shared-np  # Partition to run the job
#SBATCH --job-name=network_plot        # Job name
#SBATCH --mail-type=BEGIN              # Send email when the job starts
#SBATCH --mail-type=END                # Send email when the job ends
#SBATCH --mail-type=FAIL               # Send email if the job fails
#SBATCH --mail-user=norah.saarman@usu.edu  # Your email for notifications

# Load R module (adjust this if you need to load a specific version)
module load R  # Update with the correct R version if needed

# Run the R script
Rscript /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/scripts/network_graph.R
