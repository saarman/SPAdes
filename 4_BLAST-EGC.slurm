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
inOUTDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/spades_assembly-test/blast_result/unfiltered"
OUTDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/spades_assembly-test/blast_result/filtered_coi_out"  # Where to write output
DBDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/spades_assembly-test/blast_result/blast_DB"
COIref="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/blast_refs/coi-ref/sparrow-coi-ref.fasta"  # Reference fasta file with known sequences

# Make sure output directory exists
mkdir -p "$OUTDIR"
mkdir -p "$inOUTDIR"
mkdir -p "$DBDIR"

chmod -R g+w "$OUTDIR"
chmod -R g+w "$INDIR"
chmod -R g+w "$DBDIR"


# Loop over all filtered contig fasta files in the input directory
#for SAMPLE in `ls ./*/contigs.fasta`; do
   # Strip suffix to create a short name for output files
   # NAME=`echo $SAMPLE | sed s/_filt200-3k_sorted_contigs.fasta//g`
   # BASENAME=$(basename "$SAMPLE")  # Strip path like ./input/
   # NAME="${BASENAME%_filt200-3k_sorted_contigs.fasta}"  # Remove suffix
   # DIRNAME=$(dirname "$SAMPLE")
   # NAME=$(basename "$DIRNAME")
   # echo "Processing $NAME"
   # mkdir -p "${DBDIR}/${NAME}"

   # makeblastdb -in "$SAMPLE" -dbtype nucl
   # blastn -query "$COIref" -db ./${NAME}/contigs.fasta -out ${inOUTDIR}/${NAME}.m8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -evalue 0.05 -word_size 7 -num_threads $SLURM_CPUS_ON_NODE
   # awk '{ if ($13 > 0 && ($4 / $13) >= 0.12) print }' ${inOUTDIR}/${NAME}.m8 > ${OUTDIR}/${NAME}.m8
# done

for SAMPLE in ./*/contigs.fasta; do
   DIRNAME=$(dirname "$SAMPLE")
   NAME=$(basename "$DIRNAME")
   echo "Processing $NAME"

   # Create DB in same directory as the fasta input (required for default file name resolution)
   makeblastdb -in "$SAMPLE" -dbtype nucl -out "${DIRNAME}/contigs"

   blastn -query "$COIref" \
          -db "${DIRNAME}/contigs" \
          -out "${inOUTDIR}/${NAME}.m8" \
          -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
          -evalue 0.05 \
          -word_size 7 \
          -num_threads "$SLURM_CPUS_ON_NODE"

   awk '{ if ($13 > 0 && ($4 / $13) >= 0.12) print }' "${inOUTDIR}/${NAME}.m8" > "${OUTDIR}/${NAME}.m8"
done



chmod -R g+w "$inOUTDIR"
chmod -R g+w "$OUTDIR"
chmod -R g+w "$INDIR"
chmod -R g+w "../blast_result/"
chmod -R g+w "../SPAdes_scripts/"
