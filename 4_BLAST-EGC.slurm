#!/bin/bash
#SBATCH --time=336:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20          # same as $max set in ForkManager
#SBATCH --account=saarman-np
#SBATCH --partition=saarman-shared-np   
#SBATCH --job-name=MMseqs2_try4
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=emily.calhoun@usu.edu

# Change to BLAST directory
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/spades_assembly-test

# Load BLAST module (update if newer version becomes available)
module load blast/2.11.0

# Set input and output paths
INDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/spades_assembly-test"    # Directory with sample fasta files
inOUTDIR="$INDIR/blast_result/unfiltered"     # Raw BLAST output
OUTDIR="$INDIR/blast_result/filtered_coi_out" # Filtered output
DBDIR="$INDIR/blast_result/blast_DB"          # BLAST DB directory
COIref="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/ref/BOLD_Public.16-May-2025.fasta" #using all of BOLD as reference?
# "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/blast_refs/coi-ref/sparrow-coi-ref.fasta"  # Reference fasta file with known sequences

# Make sure directories exist
mkdir -p "$OUTDIR" "$inOUTDIR" "$DBDIR"
chmod -R g+w "$OUTDIR" "$INDIR" "$DBDIR"

# Step 1: Create BLAST database from COI reference (Only needs to be done once)

REF_DB="$DBDIR/coi_bold"

if [ ! -f "${REF_DB}.nin" ]; then
    echo "Creating BLAST database from COI reference..."
    makeblastdb -in "$COIref" -dbtype nucl -out "$REF_DB"
else
    echo "BLAST database already exists. Skipping creation."
fi

# Step 2: Loop over SPAdes contigs and BLAST them

for SAMPLE in "$INDIR"/*/contigs.fasta; do
    if [ ! -f "$SAMPLE" ]; then
        echo "No contigs.fasta files found in $INDIR. Exiting."
        exit 1
    fi

    DIRNAME=$(dirname "$SAMPLE")
    NAME=$(basename "$DIRNAME")
    echo "Processing sample: $NAME"

    # Run BLAST: contigs (query) vs COI reference DB (target)
    blastn -query "$SAMPLE" \
           -db "$REF_DB" \
           -out "${inOUTDIR}/${NAME}.m8" \
           -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
           -evalue 1e-5 \
           -word_size 5 \
           -task blastn \
           -soft_masking false \
           -dust no \
           -num_threads "${SLURM_CPUS_ON_NODE:-4}"

    # Filter hits: require alignment length ≥100 and percent identity ≥90
    awk '$4 >= 100 && $3 >= 90' "${inOUTDIR}/${NAME}.m8" > "${OUTDIR}/${NAME}.m8"
done

# Step 3: Final permissions

chmod -R g+w "$inOUTDIR" "$OUTDIR" "$INDIR"

echo "✅ BLAST search complete. Results in $OUTDIR"
