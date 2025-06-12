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
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2

# Load BLAST module (update if newer version becomes available)
module load blast/2.11.0

# Set input and output paths
INDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input"    # Directory with sample fasta files
inOUTDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/inoutputCOI"
OUTDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/outputCOI"  # Where to write output
DBDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/blastdbs"
COIref="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/ref/COIref.fasta"  # Reference fasta file with known sequences

# Make sure output directory exists
mkdir -p "$OUTDIR"
mkdir -p "$DBDIR"

chmod -R g+w "$OUTDIR"
chmod -R g+w "$INDIR"
chmod -R g+w "$DBDIR"
chmod -R g+w "../MMseqs2/"


# Loop over all filtered contig fasta files in the input directory
for SAMPLE in `ls ./input/*filt200-3k_sorted_contigs.fasta`; do
   # Strip suffix to create a short name for output files
   NAME=`echo $SAMPLE | sed s/_filt200-3k_sorted_contigs.fasta//g`
   echo "Processing $NAME"

   makeblastdb -in $SAMPLE -dbtype nucl -out ${DBDIR}/${NAME}_db
   blastn -query ./rawinput/mmRefs.fasta -db ${DBDIR}/${NAME}_db -out ${inOUTDIR}/${NAME}.m8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -max_target_seqs 10 -num_threads $(nproc)
   awk '{ if ($13 > 0 && ($4 / $13) >= 0.25) print }' ${inOUTDIR}/${NAME}.m8 > ${OUTDIR}/${NAME}.m8
done

chomd -R g+w "$inOUTDIR"
chmod -R g+w "$OUTDIR"
chmod -R g+w "$INDIR"










