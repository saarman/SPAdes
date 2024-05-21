# SPAdes
Culex species identification and blood meal analysis with de novo assembly of illumina reads from amplicon resequencing (COI, Ace2, cqm1)

# A place to keep track of Culex blood meal bioinformatics steps

## Logging onto CHPC with Terminal on a mac
1. Open Terminal
2. ssh u6036559@notchpeak.chpc.utah.edu        #replace with your username
3. salloc --time=72:00:00 --ntasks 1 --mem=100G --account=saarman-np --partition=saarman-shared-np

## Outline of steps:
1. seqtk to subsample reads (tens to hundreds of 1000's?)
2. fastqc to trim/clean/quality control
3. de novo assembly with SPAdes https://biomedicalhub.github.io/genomics/04-part4-denovo-assembly.html
4. map reads back to assembly and confirm which have good support
5. blast all consensus sequences (.fa) - relabel as COi if match COi
6. blast/bold consensus of COi matches only

## Manuals for seqtk, fastqc, SPades 
 - seqtk https://docs.csc.fi/apps/seqtk/#usage   
 - fastqc https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
 - SPAdes https://ablab.github.io/spades/ 

# Task 1: seqtk 
## working directory
```cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles```

## make dir for subseq
```mkdir subseq```

## change permissions
```chmod -R g+w subseq```

## check samples exist
```
bash
for SAMPLE in `echo B002f_S1 B013f_S2 B015f_S3 B016f_S4 B020f_S5 B021f_S6 B022f_S7 B023f_S8`; do
  echo $SAMPLE
  ls -l MAD_21_${SAMPLE}_R1_001.fastq.gz
  ls -l MAD_21_${SAMPLE}_R2_001.fastq.gz
done
```

## load the module
```module load seqtk```

## run subseq for paired reads use same seed
```
for SAMPLE in `echo B002f_S1 B013f_S2 B015f_S3 B016f_S4 B020f_S5 B021f_S6 B022f_S7 B023f_S8`; do
  echo $SAMPLE
  seqtk sample -s100 MAD_21_${SAMPLE}_R1_001.fastq.gz 10000 > ./subseq/MAD_21_${SAMPLE}_R1_001_sub1.fq; chmod -R g+w subseq
  seqtk sample -s100 MAD_21_${SAMPLE}_R2_001.fastq.gz 10000 > MAD_21_${SAMPLE}_R2_001_sub2.fq; chmod -R g+w subseq
done
```

## change permissions
```
chmod -R g+w subseq
```

# Task 2: FastQC

## example command
```
fastqc ./subseq/*.fq -o ./fastqc
```
# Task 3: SPAdes

## example command
```
module load SPAdes
spades.py -1 left.fastq.gz -2 right.fastq.gz -o output_folder
```
## make an directory   
denovo_assembly

## de novo assembly with one sample at a time
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles
spades.py -1 ./subseq/MAD_21_B002f_S1_R1_001_sub1.fq -2	./subseq/MAD_21_B002f_S1_R2_001_sub2.fq -o ./denovo_assembly
```

## de novo assembly with a loop
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles
for SAMPLE in `echo B002f_S1 B013f_S2 B015f_S3 B016f_S4 B020f_S5 B021f_S6 B022f_S7 B023f_S8`; do
  echo $SAMPLE
  spades.py -1 ./subseq/MAD_21_${SAMPLE}_R1_001_sub1.fq -2	./subseq/MAD_21_${SAMPLE}_R2_001_sub2.fq -o denovo_assembly
done
```
## change permissions
```
chmod -R g+w denovo_assembly
```




