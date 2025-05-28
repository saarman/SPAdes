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

# Change to MMSeqs2 directory
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2

# Load MMseqs2 module (update if newer version becomes available)
module load mmseqs2/oct24

# Set input and output paths
INDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input"    # Directory with sample fasta files
OUTDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/outputCOI"  # Where to write output
COIref="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/ref/COIref.fasta"  # Reference fasta file with known sequences
TEMP="/scratch/general/vast/u6036559/spades_tmp"  # Temporary directory for MMseqs2 runtime files

# Make sure output directory exists
mkdir -p "$OUTDIR"
chmod -R g+w "$OUTDIR"
chmod -R g+w "$INDIR"
chmod -R g+w "../MMseqs2/"

# Make sure the temp directory exists
mkdir -p "$TEMP"
chmod -R g+w "$TEMP"

# Loop over all filtered contig fasta files in the input directory
for SAMPLE in `ls ./input/*filt200-3k_sorted_contigs.fasta`; do
   # Strip suffix to create a short name for output files
   NAME=`echo $SAMPLE | sed s/_filt200-3k_sorted_contigs.fasta//g`
   echo "Processing $NAME"

   # Run MMseqs2 easy-search:
   # - Query = reference database (COIref)
   # - Target = sample sequence file
   # - Output = .m8 BLAST tabular format file
   # - TEMP = temporary folder for intermediate files
   # --search-type 3 = nucleotide vs nucleotide
   # --max-seqs 1 = return only the best match per query sequence
   mmseqs easy-search $coiREF $SAMPLE ${OUTDIR}/${NAME}.m8 $TEMP \
       --search-type 3 \
       --threads 20 \
       -m 8 \
done

chmod -R g+w "$OUTDIR"
chmod -R g+w "$INDIR"










