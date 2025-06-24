#!/bin/bash
#SBATCH --time=336:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20          # same as $max set in ForkManager
#SBATCH --account=saarman-np
#SBATCH --partition=saarman-shared-np   
#SBATCH --job-name=BLASTBOLD
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=emily.calhoun@usu.edu

# Load BLAST
module load blast/2.11.0

# Set paths
INDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/spades_assembly-test"
inOUTDIR="$INDIR/blast_result/unfiltered"
OUTDIR="$INDIR/blast_result/filtered_coi_out"
DBDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/spades_assembly-test/BOLD_db"
REF_DB="${DBDIR}/coi_bold"

mkdir -p "$inOUTDIR" "$OUTDIR"

# Loop through contigs and run blastn (only query vs REF DB, not creating new DBs)
for SAMPLE in "$INDIR"/*/contigs.fasta; do
    DIRNAME=$(dirname "$SAMPLE")
    NAME=$(basename "$DIRNAME")
    echo "Processing $NAME"

    blastn -query "$SAMPLE" \
           -db "$REF_DB" \
           -out "${inOUTDIR}/${NAME}.m8" \
           -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
           -evalue 10 \
           -word_size 5 \
           -task blastn \
           -soft_masking false \
           -dust no \
           -num_threads "$SLURM_CPUS_ON_NODE"

    awk '{ if ($13 > 0 && ($4 / $13) >= 0.12) print }' "${inOUTDIR}/${NAME}.m8" > "${OUTDIR}/${NAME}.m8"
done

# Set permissions (optional)
chmod -R g+w "$inOUTDIR" "$OUTDIR"
