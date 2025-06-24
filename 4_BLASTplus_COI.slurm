#!/bin/sh
#SBATCH --time=336:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20          # same as $max set in ForkManager
#SBATCH --account=saarman-np
#SBATCH --partition=saarman-shared-np   
#SBATCH --job-name=BLAST-IT!
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=A02353507@usu.edu

# Change to BLAST directory
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/denovo_assembly

# Load BLAST module (update if newer version becomes available)
module load blast/2.11.0

# Set input and output paths
INDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/denovo_assembly"    # Directory with sample fasta files
inOUTDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/blast_result/unfiltered" # Raw BLAST output
OUTDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/blast_result/filtered_coi_out"  # Filtered output
COIref="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/blast_refs/coi-ref/sparrow-coi-ref.fasta"  # Reference fasta file with sparrow sequence passer domesticus


# Make sure directories exist
mkdir -p "$OUTDIR" "$inOUTDIR"
chmod -R g+w "$OUTDIR" "$INDIR"


for SAMPLE in ./*/contigs.fasta; do
   DIRNAME=$(dirname "$SAMPLE")
   NAME=$(basename "$DIRNAME")
   echo "Processing $NAME"

   # Create DB in same directory as the fasta input (required for default file name resolution)
   makeblastdb -in "$SAMPLE" -dbtype nucl -out "${DIRNAME}/contigs"

   # Search each blast directory with Sparrow COI to identify all COIs
   blastn -query "$COIref" \
          -db "${DIRNAME}/contigs" \
          -out "${inOUTDIR}/${NAME}.m8" \
          -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
          -evalue 10 \
          -word_size 5 \
          -task blastn \
          -soft_masking false \
          -dust no \
          -num_threads "$SLURM_CPUS_ON_NODE"

   # filter blast results by coverage (0.125)
   awk '{ if ($13 > 0 && ($4 / $13) >= 0.125 && $3 > 74) print }' "${inOUTDIR}/${NAME}.m8" > "${OUTDIR}/${NAME}.m8"
done


chmod -R g+w "$inOUTDIR"
chmod -R g+w "$OUTDIR"
chmod -R g+w "$INDIR"
chmod -R g+w "../blast_result/"
chmod -R g+w "../SPAdes_scripts/"
