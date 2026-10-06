#!/bin/bash
#SBATCH --partition=saarman-np
#SBATCH --account=saarman-np
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name="readcounts"
#SBATCH --mail-user=emily.calhoun@usu.edu
#SBATCH --mail-type=BEGIN,END,FAIL

cd /uufs/chpc.utah.edu/common/home/u1055819/saarman-group/uphlfiles

find UT-* -name "*_R1_001.fastq.gz" | grep -v sub | \
  xargs -P 8 -I{} bash -c \
  'n=$(( $(zcat "{}" | wc -l) / 4 )); echo -e "$(dirname "{}")\t$(basename "{}" _R1_001.fastq.gz)\t$n"' \
  > ~/read_counts_full.tsv
