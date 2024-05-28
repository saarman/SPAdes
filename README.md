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
3. SPAdes for de novo assembly, https://biomedicalhub.github.io/genomics/04-part4-denovo-assembly.html
4. Filter contigs for min length and min coverage
5. Blastn to identify COi match

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
bash
for SAMPLE in `echo B002f_S1 B013f_S2 B015f_S3 B016f_S4 B020f_S5 B021f_S6 B022f_S7 B023f_S8`; do
  echo $SAMPLE
  spades.py -1 ./subseq/MAD_21_${SAMPLE}_R1_001_sub1.fq -2	./subseq/MAD_21_${SAMPLE}_R2_001_sub2.fq -o ./denovo_assembly/${SAMPLE} --isolate
done
```
## change permissions
```
chmod -R g+w denovo_assembly
```

# Step 3: Filter for minimum coverage and length

## With a hard threshold filter of COVERAGE>600, LENGTH> 200 
Remove unreliable samples (no COi amplification, according to tape station results). 
B0021-23 did not have successful COi PCR amplification. 
B015f did have some amplification with only 602 coverage, so will place threshold at 600. 
B023f did not work, but had maximum 531 coverage on a random fragment.
This threshold depends on the subsampling step as well, so needs to be adjusted for each experiment. 

```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/denovo_assembly
bash
COVERAGE=600 #change this to set min coverage
LENGTH=200 #change this to set min length
for SAMPLE in `echo B002f_S1 B013f_S2 B015f_S3 B016f_S4 B020f_S5 B021f_S6 B022f_S7 B023f_S8`; do
  echo $SAMPLE
  perl -0076 -ne 'chomp;($h,@S)=split/\n/;$s=join("",@S);print"$h\t$s\n"unless(!$h)' ./${SAMPLE}/contigs.fasta | sed 's/_/ /g' | awk -F " " -v a="$LENGTH" -v b="$COVERAGE" '$4>=a && $6>=b' | sed 's/ /_/g' | sed 's/\t/\n/g' | sed "s/NODE/\>${SAMPLE}/g" > ./${SAMPLE}/filtered600_contigs.fasta
 cp ./${SAMPLE}/filtered600_contigs.fasta ./${SAMPLE}_filtered_contigs.fasta
done
chmod -R g+w ../denovo_assembly
cat *.fasta #concatenate
```
* I also noticed that some of the sequences had poly-C and poly-G at ends, probably an artifact that can be removed by filtering raw reads before assembly.

## With a SORT for top 3 COVERAGE added 
Still includes hard filter for LENGTH>200 and COVERAGE>30

```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/denovo_assembly
bash
LENGTH=200 #change this to set min length
COVERAGE=30 #change this to set min coverage
for SAMPLE in `echo B002f_S1 B013f_S2 B015f_S3 B016f_S4 B020f_S5 B021f_S6 B022f_S7 B023f_S8`; do
  echo $SAMPLE
  perl -0076 -ne 'chomp;($h,@S)=split/\n/;$s=join("",@S);print"$h\t$s\n"unless(!$h)' ./${SAMPLE}/contigs.fasta  | sed 's/_/ /g' | awk -F " " -v a="$LENGTH" -v b="$COVERAGE" '$4>=a && $6>=b' |  sed 's/ /_/g'  | sort -r -t_ -nk6 | head -3 | sed 's/\t/\n/g' | sed "s/NODE/\>${SAMPLE}/g" > ./${SAMPLE}/sort30_contigs.fasta
 cp ./${SAMPLE}/sort30_contigs.fasta ./${SAMPLE}_sorted_contigs.fasta
done
chmod -R g+w ../denovo_assembly
cat *sorted_contigs.fasta    #concatenate all sorted matches
```
sort -nk3

# Step 4: Blastn
