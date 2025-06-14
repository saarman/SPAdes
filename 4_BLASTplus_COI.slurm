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

# Change to BLAST directory
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/denovo_assembly

# Load BLAST module (update if newer version becomes available)
module load blast/2.11.0

# Set input and output paths
INDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/denovo_assembly"    # Directory with sample fasta files
inOUTDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/ej-spades-scripts/blast_result/unfiltered"
OUTDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/ej-spades-scripts/blast_result/filtered_coi_out"  # Where to write output
DBDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/ej-spades-scripts/blast_result/blast_DB"
COIref="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/blast_refs/coi-ref/sparrow-coi-ref.fasta"  # Reference fasta file with known sequences

# Make sure output directory exists
mkdir -p "$OUTDIR"
mkdir -p "$inOUTDIR"
mkdir -p "$DBDIR"

chmod -R g+w "$OUTDIR"
chmod -R g+w "$INDIR"
chmod -R g+w "$DBDIR"
chmod -R g+w "../denovo_assembly/"
chmod -R g+w "../blast_result/"

# Loop over all filtered contig fasta files in the input directory
for SAMPLE in `ls ./*/contigs.fasta`; do
   # Strip suffix to create a short name for output files
   # NAME=`echo $SAMPLE | sed s/_filt200-3k_sorted_contigs.fasta//g`
   # BASENAME=$(basename "$SAMPLE")  # Strip path like ./input/
   # NAME="${BASENAME%_filt200-3k_sorted_contigs.fasta}"  # Remove suffix
   DIRNAME=$(dirname "$SAMPLE")
   NAME=$(basename "$DIRNAME")
   echo "Processing $NAME"
   mkdir -p "${DBDIR}/${NAME}"

   makeblastdb -in "$SAMPLE" -dbtype nucl -out ${DBDIR}/${NAME}/${NAME}_db
   blastn -query "$COIref" -db ${DBDIR}/${NAME}/${NAME}_db -out ${inOUTDIR}/${NAME}.m8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -max_target_seqs 10 -word_size 7 -num_threads $SLURM_CPUS_ON_NODE
   awk '{ if ($13 > 0 && ($4 / $13) >= 0.25) print }' ${inOUTDIR}/${NAME}.m8 > ${OUTDIR}/${NAME}.m8
done

chmod -R g+w "$inOUTDIR"
chmod -R g+w "$OUTDIR"
chmod -R g+w "$INDIR"










