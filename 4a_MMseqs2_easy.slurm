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

# Load modules
module load mmseqs2/oct24  # change to module name

# Assign variables – inputs and outputs
bash
INDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input"
OUTDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output"
REF="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input/mmREF.fasta"
TEMP="/scratch/general/vast/u6036559/spades_tmp/"

# Run Command in a loop:
cd $INDIR
for SAMPLE in `ls *filt200-3k_sorted_contigs.fasta`; do
   NAME=`echo $SAMPLE | sed s/_filt200-3k_sorted_contigs.fasta//g`
   echo $NAME
   mmseqs easy-search $REF $SAMPLE ${OUTDIR}/${NAME}.m8 $TEMP --search-type 3 --threads 20
done
