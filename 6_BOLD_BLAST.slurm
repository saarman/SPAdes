#!/bin/bash
#SBATCH --time=336:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20          # same as $max set in ForkManager
#SBATCH --account=saarman-np
#SBATCH --partition=saarman-shared-np   
#SBATCH --job-name=consensus-COI
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=emily.calhoun@usu.edu


# Change to BLAST directory
# cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/denovo_assembly

# Load any needed modules 
module load blast/2.11.0


# Set input and output paths
INFILE="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/blast_result/consensus_multi_hits.fasta"    # Directory with sample fasta files
DataBASE="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/BOLD_Chordata_db/BOLDChordata"
OUTDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/blast_result"  # Filtered output


# Make sure directories exist
mkdir -p "$OUTDIR"
chmod -R g+w "$OUTDIR"


# Search each blast directory with Sparrow COI to identify all COIs
   blastn -query "$INFILE" \
          -db "$DataBASE" \
          -out "${OUTDIR}/test-bold-chordata-hits.m8" \
          -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
          -evalue 1e-10 \
          -perc_identity 95 \
          -word_size 5 \
          -max_target_seqs 1000000 \
          -task blastn \
          -soft_masking false \
          -dust no \
          -num_threads "$SLURM_CPUS_ON_NODE"




# awk '{ if ($13 > 0 && ($4 / $13) >= 0.125 && $3 > 74) print }' "${inOUTDIR}/${NAME}.m8" > "${OUTDIR}/${NAME}.m8"



