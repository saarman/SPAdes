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
#SBATCH --mail-user=emily.calhoun@usu.edu

# Change to BLAST directory
# cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/denovo_assembly

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
chmod -R g+w "../SPAdes_scripts/"
chmod -R g+w "../blast_result/"

# Loop over all filtered contig fasta files in the input directory
for SAMPLE in `ls ./input/*filt200-3k_sorted_contigs.fasta`; do
   # Strip suffix to create a short name for output files
   # NAME=`echo $SAMPLE | sed s/_filt200-3k_sorted_contigs.fasta//g`
   BASENAME=$(basename "$SAMPLE")  # Strip path like ./input/
   NAME="${BASENAME%_filt200-3k_sorted_contigs.fasta}"  # Remove suffix
   echo "Processing $NAME"
   mkdir -p "${DBDIR}/${NAME}"

   makeblastdb -in $SAMPLE -dbtype nucl -out ${DBDIR}/${NAME}/${NAME}_db
   blastn -query ./rawinput/mmRefs.fasta -db ${DBDIR}/${NAME}_db -out ${inOUTDIR}/${NAME}.m8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -max_target_seqs 10 -num_threads $SLURM_CPUS_ON_NODE
   awk '{ if ($13 > 0 && ($4 / $13) >= 0.25) print }' ${inOUTDIR}/${NAME}.m8 > ${OUTDIR}/${NAME}.m8
done

chmod -R g+w "$inOUTDIR"
chmod -R g+w "$OUTDIR"
chmod -R g+w "$INDIR"
