#!/bin/bash
#SBATCH --time=336:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20          # same as $max set in ForkManager
#SBATCH --account=saarman-np
#SBATCH --partition=saarman-shared-np   
#SBATCH --job-name=BLAST-IT!
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=emily.calhoun@usu.edu



# take only top hit for each sample
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/blast_result/filtered_coi_out
for SEQ in `cat bloodmeal-coi-HDAWWRMX016-HitTable.csv | awk -F"," '{print $1}' | uniq`; do
   grep ${SEQ} -m 1 bloodmeal-coi-HDAWWRMX016-HitTable.csv >> bloodmeal-coi-HDAWWRMX016-HitTable-top1.csv
done

# write top hit and scientific name to file 
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/blast_result/filtered_coi_out
for SEQ in `cat bloodmeal-coi-HDAWWRMX016-HitTable.csv | awk -F"," '{print $1}' | uniq`; do
   HIT=`grep ${SEQ} -A 5 HDAWWRMX016-Alignment.txt | tail -1 | awk '{print $1 " " $2}'`
   TAB=`grep ${SEQ} -m 1 bloodmeal-coi-HDAWWRMX016-HitTable.csv`
   echo ${TAB}, $HIT >> bloodmeal-coi-HDAWWRMX016-HitTable-top1-readable.csv
done
