# SPAdes
Culex species identification and blood meal analysis with de novo assembly of illumina reads from amplicon resequencing (COI, Ace2, cqm1)

# A place to keep track of Culex blood meal bioinformatics steps

## Logging onto CHPC with Terminal on a mac
1. Open Terminal
2. ssh u6036559@notchpeak.chpc.utah.edu  #replace with your username
3. salloc --time=72:00:00 --ntasks 1 --mem=100G --account=saarman-np --partition=saarman-shared-np

## Outline of steps:
1. fastqc to trim/clean/quality control
2. SPAdes for de novo assembly, https://biomedicalhub.github.io/genomics/04-part4-denovo-assembly.html
3. Filter contigs for min length and min coverage
4. Clustering/Blastn to identify ace2/cqm1 haplotypes and COi match

## Manuals for fastqc, SPades 
 - fastqc https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
 - SPAdes https://ablab.github.io/spades/ 

# Task 1: FastQC

## Make output directory
```
# CD to working directory
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M70330-240718

# make dir for fastqc
mkdir fastqc

#change permissions
chmod -R g+w fastqc
```

## Run command
```
# load module 
module load fastqc

# run fastqc for all .fastq.gz files in directory
fastqc *.fastq.gz -o ./fastqc

#change permissions
chmod -R g+w fastqc
```

# Task 2: SPAdes

## Example command
module load spades
spades.py -1 left.fastq.gz -2 right.fastq.gz -o output_folder

## Make an directory and change permissions
```
# change directory
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M70330-240718

# make new directory named 'denovo_assembly'
mkdir denovo_assembly

# change permissions
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M70330-240718
```

## Run SPAdes de novo assembly with a loop

List of samples to run:
```
ls -l *R1_001.fastq.gz | awk '{print $NF}' | cut -d_ -f1-2
```

SPAdes:
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

# Step 3: Filter and sort contigs output per sample, add sample to fasta header

## Example with filters (length, coverage), sorted by coverage
Filters, sorts, and adds name of sample to header of each contig
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702/denovo_assembly
bash
LENGTH=450 #change this to set min length
COVERAGE=20 #change this to set min coverage
MAXLEN=700 #change this to set max length
for SAMPLE in `ls -l | grep -v "total" |  grep -v "fasta" | awk '{print $NF}'`; do
  echo $SAMPLE
  perl -0076 -ne 'chomp;($h,@S)=split/\n/;$s=join("",@S);print"$h\t$s\n"unless(!$h)' ./${SAMPLE}/contigs.fasta  | sed 's/_/ /g' | awk -F " " -v a="$LENGTH" -v b="$COVERAGE" -v c="$MAXLEN" '$4>=a && $6>=b && $4<=c' |  sed 's/ /_/g'  | sort -r -t_ -nk6 | sed 's/\t/\n/g' | sed "s/NODE/\>${SAMPLE}/g" > ./${SAMPLE}/contigs.fasta
 cp ./${SAMPLE}/sorted_contigs.fasta ./${SAMPLE}_sorted_contigs.fasta
done

# update permissions
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M70330-240718

# list all sorted matches
ls -l *sorted_contigs.fasta

# concatenate all sorted matches
cat *sorted_contigs.fasta    
```

# Step 4: 

