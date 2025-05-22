# SPAdes
Culex species identification and blood meal analysis with de novo assembly of illumina reads from amplicon resequencing (COI, Ace2, cqm1)

# A place to keep track of Culex blood meal bioinformatics steps

## Logging onto CHPC with Terminal on a mac
1. Open Terminal
2. ssh u6036559@notchpeak.chpc.utah.edu  #replace with your username
3. salloc --time=336:00:00 --ntasks 1 --mem=100G --account=saarman-np --partition=saarman-shared-np

## Logging onto CHPC with command line on a Windows
1. Go to Interactive Desktop at https://ondemand.chpc.utah.edu/pun/sys/dashboard/batch_connect/sys/desktop_expert/session_contexts/new
2. Cluster: notchpeak, Account and partition: saarman-np:saarman-shared-np, number of cores: 1, number of hours: 72-336, memory per job in GB: 128

## Outline of steps:
1. fastqc to trim/clean/quality control
2. SPAdes for de novo assembly, https://biomedicalhub.github.io/genomics/04-part4-denovo-assembly.html
3. Filter and sort contigs for min length and min coverage
4. Clustering/Blastn to identify ace2/cqm1 haplotypes and COi match
5. Final QC (mapping reads, blastn, determine haplotypes, etc.)
6. Blastn or BOLD

## Manuals for fastqc, SPades 
 - fastqc https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
 - SPAdes https://ablab.github.io/spades/ 

# Task 1: FastQC

Make output directory
```
# CD to working directory
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M70330-240718

# make dir for fastqc
mkdir fastqc

#change permissions
chmod -R g+w fastqc
```

Run command
```
# load module 
module load fastqc

# run fastqc for all .fastq.gz files in directory
fastqc *.fastq.gz -o ./fastqc

#change permissions
chmod -R g+w fastqc
```

# Task 2: SPAdes

Make an output directory:
```
# change directory
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M70330-240718

# make new directory named 'denovo_assembly'
mkdir denovo_assembly

# change permissions
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M70330-240718
```

List of samples to run:
```
ls -l *R1_001.fastq.gz | awk '{print $NF}' | cut -d_ -f1-2
```

Help menu:
```
module load spades
spades.py --help
```

SPAdes in a loop:
```
# change directory
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M70330-240718

# load module
module load spades

# loop
bash
for SAMPLE in `ls -l *R1_001.fastq.gz | awk '{print $NF}' | cut -d_ -f1-2`; do
  echo $SAMPLE
  spades.py -1 ${SAMPLE}_L001_R1_001.fastq.gz -2	${SAMPLE}_L001_R2_001.fastq.gz -o ./denovo_assembly/${SAMPLE} --isolate
done
```

Remember to change permissions:
```
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M70330-240718
```

# Step 3: Filter and sort contigs

Example with filters (length, coverage), sorted by coverage, adds sample ID to header of each contig
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M70330-240718/denovo_assembly
bash
LENGTH=200 #change this to set min length
COVERAGE=2 #change this to set min coverage
MAXLEN=3000 #change this to set max length
for SAMPLE in `ls -l | grep -v "total" |  grep -v "fasta" | awk '{print $NF}'`; do
  echo $SAMPLE
  perl -0076 -ne 'chomp;($h,@S)=split/\n/;$s=join("",@S);print"$h\t$s\n"unless(!$h)' ./${SAMPLE}/contigs.fasta | awk -F "_" -v a="$LENGTH" -v b="$COVERAGE" -v c="$MAXLEN" '$4>=a && $6>=b && $4<=c' | sort -r -t_ -nk6 | sed 's/\t/\n/g' | sed "s/NODE/\>${SAMPLE}/g" > ./${SAMPLE}/${SAMPLE}_filt200-3k_sorted_contigs.fasta
  cp ./${SAMPLE}/${SAMPLE}_filt200-3k_sorted_contigs.fasta /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input/${SAMPLE}_filt200-3k_sorted_contigs.fasta
done

cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702/denovo_assembly
bash
LENGTH=200 #change this to set min length
COVERAGE=2 #change this to set min coverage
MAXLEN=3000 #change this to set max length
for SAMPLE in `ls -l | grep -v "total" |  grep -v "fasta" | awk '{print $NF}'`; do
  echo $SAMPLE
  perl -0076 -ne 'chomp;($h,@S)=split/\n/;$s=join("",@S);print"$h\t$s\n"unless(!$h)' ./${SAMPLE}/contigs.fasta | awk -F "_" -v a="$LENGTH" -v b="$COVERAGE" -v c="$MAXLEN" '$4>=a && $6>=b && $4<=c' | sort -r -t_ -nk6 | sed 's/\t/\n/g' | sed "s/NODE/\>${SAMPLE}/g" > ./${SAMPLE}/${SAMPLE}_filt200-3k_sorted_contigs.fasta
  cp ./${SAMPLE}/${SAMPLE}_filt200-3k_sorted_contigs.fasta /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input/${SAMPLE}_filt200-3k_sorted_contigs.fasta
done

# update permissions
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/

# list all sorted matches
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input/; ls -l *.fasta

# concatenate all sorted matches
cat /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input/*.fasta | head -20000
```

# Step 4: Search and Clustering with MMseqs2

## Step 4a Easy-Search with mmRef.fasta vs each .fasta
Example of easy-search commands in 4a_MMseqs2_easy.slurm
```
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
REF="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/scripts/mmRefs.fasta"
TEMP="/scratch/general/vast/u6036559/spades_tmp/"

# Run Command in a loop:
cd $INDIR
for SAMPLE in `ls *filt200-3k_sorted_contigs.fasta`; do
   NAME=`echo $SAMPLE | sed s/_filt200-3k_sorted_contigs.fasta//g`
   echo $NAME
   mmseqs easy-search $REF $SAMPLE ${OUTDIR}/${NAME}.m8 $TEMP --search-type 3 --threads 20
done
```


## Step 4b Pull out best e-val for each reference for each .fasta 
  
EF061759.1 cqm1   
AY497524.1 ace2    
AY666266.1 + GU130589.1 COi (house finch, cow)    

Example of extract commands in step 4b:
```
salloc --time=336:00:00 --ntasks 1 --mem=100G --account=saarman-np --partition=saarman-shared-np

# Assign variables – inputs and outputs
bash
INDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input"
OUTDIR="/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output"

# Run Command in a loop for cqm1:
cd $OUTDIR
for SAMPLE in `ls B*.m8`; do
   NAME=`echo $SAMPLE | sed s/.m8//g`
   echo $NAME
   SEQ=`cat ${NAME}.m8 | grep -m 1 "EF061759.1" | awk '{print $2}'`
   cat ${INDIR}/${NAME}*.fasta | grep $SEQ -A 1 >> cqm1_matches.fasta
done

# Run Command in a loop for ace2:
cd $OUTDIR
for SAMPLE in `ls B*.m8`; do
   NAME=`echo $SAMPLE | sed s/.m8//g`
   echo $NAME
   SEQ=`cat ${NAME}.m8 | grep -m 1 "AY497524.1" | awk '{print $2}'`
   cat ${INDIR}/${NAME}*.fasta | grep $SEQ -A 1 >> ace2_matches.fasta
done

# Run Command in a loop for COi:
cd $OUTDIR
for SAMPLE in `ls B*.m8`; do
   NAME=`echo $SAMPLE | sed s/.m8//g`
   echo $NAME
   SEQ1=`cat ${NAME}.m8 | grep -m 1 "AY666266.1"| awk '{print $2}'`
   SEQ2=`cat ${NAME}.m8 | grep -m 1 "GU130589.1"| awk '{print $2}'`
   if [[ $SEQ1 == $SEQ2 ]]; then
      cat ${INDIR}/${NAME}*.fasta | grep $SEQ1 -A 1 >> coi_matches.fasta
   else
      cat ${INDIR}/${NAME}*.fasta | grep $SEQ1 -A 1 >> coi_matches.fasta
      cat ${INDIR}/${NAME}*.fasta | grep $SEQ2 -A 1 >> coi_matches.fasta
   fi
done


```

## Github
***Run just ONCE: Clone from Github***  
Before running, I need to make these files on github, and then use git to clone
```
# Just once:
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2
git clone https://github.com/saarman/SPAdes scripts
```
***Need to run every time: Pull from Github and run with Sbatch***  
Pull and run with sbatch
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/scripts
git pull
sbatch 4a_MMseqs2_easy.slurm
```

# Step 5. QC and map reads to de novo contigs
Placeholder for now

# Step 6. Blastn or BOLD Systems
NOTE: I tried using the BOLD Systems search engine, but this accepts only forward strand, you would need to reverse-complement a random assortment before searching... not convenient.

First I ran online ncbi blastn with the coi_matches.fasta file created above in **Step 4b**   

Then I placed results on chpc in these locations:   
/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/blastn/bloodmeal-coi-HDAWWRMX016-HitTable.csv  
/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/blastn/HDAWWRMX016-Alignment.txt  
FORMAT: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6, as well as the alignments that give the scientific name of the top hit.  

## Filter blastn results for top hit for each sequence plus scientific name
```
# update scripts
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/scripts
git pull

# change permissions
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/blastn/

# take only top hit for each sample
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/blastn/
for SEQ in `cat bloodmeal-coi-HDAWWRMX016-HitTable.csv | awk -F"," '{print $1}' | uniq`; do
   grep ${SEQ} -m 1 bloodmeal-coi-HDAWWRMX016-HitTable.csv >> bloodmeal-coi-HDAWWRMX016-HitTable-top1.csv
done

# write top hit and scientific name to file 
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/blastn/
for SEQ in `cat bloodmeal-coi-HDAWWRMX016-HitTable.csv | awk -F"," '{print $1}' | uniq`; do
   HIT=`grep ${SEQ} -A 5 HDAWWRMX016-Alignment.txt | tail -1 | awk '{print $1 " " $2}'`
   TAB=`grep ${SEQ} -m 1 bloodmeal-coi-HDAWWRMX016-HitTable.csv`
   echo ${TAB}, $HIT >> bloodmeal-coi-HDAWWRMX016-HitTable-top1-readable.csv
done
```
  
# Step 7: Deal with  multiple feedings... unclear what is the best method 

Top several hits?  
Threshold percent identity?   
Run step 4b again with birds/mammal ref separately first? Then take top hit(s) and apply threshold?   







