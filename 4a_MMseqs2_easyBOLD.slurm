#!/bin/sh
#SBATCH --time=336:00:00                  # Maximum run time (14 days)
#SBATCH --nodes=1                         # Run on a single node
#SBATCH --ntasks=20                       # Use 20 CPU threads
#SBATCH --account=saarman-np              # CHPC account
#SBATCH --partition=saarman-shared-np     # Partition/queue to submit job
#SBATCH --job-name=MMseqs2_try4           # Job name for queue tracking
#SBATCH --mail-type=BEGIN,END,FAIL        # Email notifications for job events
#SBATCH --mail-user=norah.saarman@usu.edu # Your email address for notifications

# Load MMseqs2 module (update if newer version becomes available)
module load mmseqs2/oct24

# Set input and output paths
INDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input"    # Directory with sample fasta files
OUTDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/outputBOLD"  # Where to write output
REF="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/ref/BOLD_Public.16-May-2025.fasta"  # Reference fasta file with known sequences
TEMP="/scratch/general/vast/u6036559/spades_tmp/"  # Temporary directory for MMseqs2 runtime files

# Make sure output directory exists
mkdir -p "$OUTDIR"
chmod -R g+w "$OUTDIR"
chmod -R g+w "$INDIR"
chmod -R g+w "../MMseqs2/"

# Make sure the temp directory exists
TEMP="${OUTDIR}/mmseqs_tmp_${NAME}"
mkdir -p "$TEMP"
chmod -R g+w "$TEMP"

# Change to input directory
cd "$INDIR"

# Loop over all filtered contig fasta files in the input directory
for SAMPLE in `ls *filt200-3k_sorted_contigs.fasta`; do
   # Strip suffix to create a short name for output files
   NAME=`echo $SAMPLE | sed s/_filt200-3k_sorted_contigs.fasta//g`
   echo "Processing $NAME"

   # Run MMseqs2 easy-search:
   # - Query = reference database (REF)
   # - Target = sample sequence file
   # - Output = .m8 BLAST tabular format file
   # - TEMP = temporary folder for intermediate files
   # --search-type 3 = nucleotide vs nucleotide
   # --max-seqs 1 = return only the best match per query sequence
   mmseqs easy-search $REF $SAMPLE ${OUTDIR}/${NAME}.m8 $TEMP \
       --search-type 3 \
       --threads 20 \
       --max-seqs 1
done

chmod -R g+w "$OUTDIR"
chmod -R g+w "$INDIR"
