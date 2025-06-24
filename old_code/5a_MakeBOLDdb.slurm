#!/bin/bash
#SBATCH --time=336:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20          # same as $max set in ForkManager
#SBATCH --account=saarman-np
#SBATCH --partition=saarman-shared-np   
#SBATCH --job-name=MakeBOLDdb
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=emily.calhoun@usu.edu

module load blast/2.11.0

COIref="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/BOLD_Chordata_db"
DBDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/BOLD_Chordata_db"
REF_DB="${DBDIR}/coi_bold"

mkdir -p "$DBDIR"

if [ ! -f "${REF_DB}.nin" ]; then
    echo "Creating BLAST database from BOLD reference..."
    makeblastdb -in "$COIref" -dbtype nucl -out "$REF_DB"
else
    echo "BLAST DB already exists. Skipping."
fi
