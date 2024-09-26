# July 2024
## new data from UPHL, first P1 2018 and P1 2021
salloc --time=72:00:00 --ntasks 1 --mem=100G --account=saarman-np --partition=saarman-shared-np

## Run md5sum
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702
md5sum *fastq.gz > received_checksum.txt

## Check for differences before and after transfer
diff md5sum.txt  received_checksum.txt 

## Make list of names of samples we want to loop through:
ls -l *R1_001.fastq.gz |  awk '{print $NF}' | cut -d_ -f1-2

## Apply this to a for loop:
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702
for SAMPLE in `ls -l *R1_001.fastq.gz |  awk '{print $NF}' | cut -d_ -f1-2`; do
  echo $SAMPLE
  spades.py -1 ${SAMPLE}_L001_R1_001.fastq.gz -2 ${SAMPLE}_L001_R2_001.fastq.gz -o ./denovo_assembly/${SAMPLE} --isolate
done


# 09/26/2024 Dump of pipeline so far: 

# SPAdes
Culex species identification and blood meal analysis with de novo assembly of illumina reads from amplicon resequencing (COI, Ace2, cqm1)

## Outline of steps:
1. (Optional?) seqtk to subsample reads (tens to hundreds of 1000's?)
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

## Trial Step 4: Filter for minimum coverage and length

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
NUM=3 #change this to set number of top sorted contigs to retain
for SAMPLE in `echo B002f_S1 B013f_S2 B015f_S3 B016f_S4 B020f_S5 B021f_S6 B022f_S7 B023f_S8`; do
  echo $SAMPLE
  perl -0076 -ne 'chomp;($h,@S)=split/\n/;$s=join("",@S);print"$h\t$s\n"unless(!$h)' ./${SAMPLE}/contigs.fasta  | sed 's/_/ /g' | awk -F " " -v a="$LENGTH" -v b="$COVERAGE" '$4>=a && $6>=b' |  sed 's/ /_/g'  | sort -r -t_ -nk6 | head -${NUM} | sed 's/\t/\n/g' | sed "s/NODE/\>${SAMPLE}/g" > ./${SAMPLE}/sort30_contigs.fasta
 cp ./${SAMPLE}/sort30_contigs.fasta ./${SAMPLE}_sorted_contigs.fasta
done
chmod -R g+w ../denovo_assembly
cat *sorted_contigs.fasta    #concatenate all sorted matches
```
sort -nk3

## Another example with SORT min coverage 150
```
bash
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702/denovo_assembly
LENGTH=150 #change this to set min length
COVERAGE=1 #change this to set min coverage
TOP=500 #change this to set number of top sorted contigs to retain
for SAMPLE in `ls -l | grep -v "total" |  grep -v "fasta" | awk '{print $NF}'`; do
  echo $SAMPLE
  cat ./${SAMPLE}/contigs.fasta | perl -0076 -ne 'chomp;($h,@S)=split/\n/;$s=join("",@S);print"$h\t$s\n"unless(!$h)' | sed 's/_/ /g' | awk -F " " -v a="$LENGTH" -v b="$COVERAGE" '$4>=a && $6>=b' |  sed 's/ /_/g'  | sort -r -t_ -nk6 | sed 's/\t/\n/g' | sed "s/NODE/\>${SAMPLE}_NODE/g" > ./${SAMPLE}_min${LENGTH}_sorted.fasta
done

# return all 
cat *sorted.fasta
#this returns 179,792 sequences to blast. Too many?
#yes, Your query(45,700,645 bases) is longer than the maximum allowed(1,000,000 bases).

#write to file to download and upload to blastn
cat *sorted.fasta > all_min150_sorted.fasta 
```

## Another example with SORT for top 10 
Still includes hard filters

```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702/denovo_assembly
bash

LENGTH=450 #change this to set min length
COVERAGE=20 #change this to set min coverage
MAXLEN=700 #change this to set max length
NUM=10 #change this to set number of top sorted contigs to retain
for SAMPLE in `ls -l | grep -v "total" |  grep -v "fasta" | awk '{print $NF}'`; do
  echo $SAMPLE
  perl -0076 -ne 'chomp;($h,@S)=split/\n/;$s=join("",@S);print"$h\t$s\n"unless(!$h)' ./${SAMPLE}/contigs.fasta  | sed 's/_/ /g' | awk -F " " -v a="$LENGTH" -v b="$COVERAGE" -v c="$MAXLEN" '$4>=a && $6>=b && $4<=c' |  sed 's/ /_/g'  | sort -r -t_ -nk6 | head -${NUM} | sed 's/\t/\n/g' | sed "s/NODE/\>${SAMPLE}/g" > ./${SAMPLE}/top${NUM}_contigs.fasta
 cp ./${SAMPLE}/top${NUM}_contigs.fasta ./${SAMPLE}_top${NUM}_contigs.fasta
done
chmod -R g+w ../denovo_assembly
ls -l *top${NUM}_contigs.fasta  #list all sorted matches
cat *top10_contigs.fasta    #concatenate all sorted matches

```
# Step 4: A different approach - filtering with primer sequence  
Idea and example at https://bioinf.shenwei.me/seqkit/usage/#seqkit  
Also included a sorting and filter step
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702/denovo_assembly
mkdir ../seqkit
chmod -R g+w ../seqkit  
```
## For COI:  
Primers for coi: Mod_RepCOI_F TNTTYTCMACYAACCACAAAGA, Mod_RepCOI_R TTCDGGRTGNCCRAARAATCA, 	VertCOI_7216_R CARAAGCTYATGTTRTTYATDCG, VertCOI_7194_F CGMATRAAYAAYATRAGCTTCTGAY   
Search for primer sequence after trimming 2+ bp from start/end, resolving as many degenerate bases as possible, with 2? mismatches allowed

*NOTE: [ERRO] flag -r (--use-regexp) or (--degenerate) not allowed when giving flag -m (--max-mismatch)*
*I found that using -m worked better anyway!*


Mod_RepCOI_F, removed 2 from start, no degenerate sites
TTCTCAACCAACCACAAAGA
TTCTCAACTAACCACAAAGA
TTCTCCACCAACCACAAAGA
TTCTCCACTAACCACAAAGA
TTTTCAACCAACCACAAAGA
TTTTCAACTAACCACAAAGA
TTTTCCACCAACCACAAAGA
TTTTCCACTAACCACAAAGA

Mod_RepCOI_R, removed 4 from start, one degenerate site
GGATGNCCAAAAAATCA
GGGTGNCCAAAAAATCA
GGATGNCCGAAAAATCA
GGGTGNCCGAAAAATCA
GGATGNCCAAAGAATCA
GGGTGNCCAAAGAATCA
GGATGNCCGAAGAATCA
GGGTGNCCGAAGAATCA

VertCOI_7216_R, removed 3 from end, no degenerate sites
CAAAAGCTCATGTTATTCAT
CAAAAGCTCATGTTATTTAT
CAAAAGCTCATGTTGTTCAT
CAAAAGCTCATGTTGTTTAT
CAGAAGCTCATGTTATTCAT
CAGAAGCTCATGTTATTTAT
CAGAAGCTCATGTTGTTCAT
CAGAAGCTCATGTTGTTTAT
CAAAAGCTTATGTTATTCAT
CAAAAGCTTATGTTATTTAT
CAAAAGCTTATGTTGTTCAT
CAAAAGCTTATGTTGTTTAT
CAGAAGCTTATGTTATTCAT
CAGAAGCTTATGTTATTTAT
CAGAAGCTTATGTTGTTCAT
CAGAAGCTTATGTTGTTTAT

VertCOI_7194_F, removed 3 from start, 1 from end, no degenerate sites
ATAAACAACATAAGCTTCTGA
ATAAATAACATAAGCTTCTGA
ATAAACAACATGAGCTTCTGA
ATAAATAACATGAGCTTCTGA
ATAAACAATATAAGCTTCTGA
ATAAATAATATAAGCTTCTGA
ATAAACAATATGAGCTTCTGA
ATAAATAATATGAGCTTCTGA
ATGAACAACATAAGCTTCTGA
ATGAATAACATAAGCTTCTGA
ATGAACAACATGAGCTTCTGA
ATGAATAACATGAGCTTCTGA
ATGAACAATATAAGCTTCTGA
ATGAATAATATAAGCTTCTGA
ATGAACAATATGAGCTTCTGA
ATGAATAATATGAGCTTCTGA



```
bash
module load seqkit
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702/denovo_assembly
mkdir ../seqkit/coi ; chmod -R g+w ../seqkit/coi
LENGTH=150 #change this to set min length
COVERAGE=1 #change this to set min coverage
TOP=1 #change this to set number of top sorted contigs to retain

for SAMPLE in `ls -l | grep -v "total" |  grep -v "fasta" | awk '{print $NF}'`; do
  echo $SAMPLE
   cat ./${SAMPLE}/contigs.fasta | seqkit grep -s -i -p TTCTCAACCAACCACAAAGA -p TTCTCAACTAACCACAAAGA -p TTCTCCACCAACCACAAAGA -p TTCTCCACTAACCACAAAGA -p TTTTCAACCAACCACAAAGA -p TTTTCAACTAACCACAAAGA -p TTTTCCACCAACCACAAAGA -p TTTTCCACTAACCACAAAGA -p GGATGNCCAAAAAATCA -p GGGTGNCCAAAAAATCA -p GGATGNCCGAAAAATCA -p GGGTGNCCGAAAAATCA -p GGATGNCCAAAGAATCA -p GGGTGNCCAAAGAATCA -p GGATGNCCGAAGAATCA -p GGGTGNCCGAAGAATCA -p CAAAAGCTCATGTTATTCAT -p CAAAAGCTCATGTTATTTAT -p CAAAAGCTCATGTTGTTCAT -p CAAAAGCTCATGTTGTTTAT -p CAGAAGCTCATGTTATTCAT -p CAGAAGCTCATGTTATTTAT -p CAGAAGCTCATGTTGTTCAT -p CAGAAGCTCATGTTGTTTAT -p CAAAAGCTTATGTTATTCAT -p CAAAAGCTTATGTTATTTAT -p CAAAAGCTTATGTTGTTCAT -p CAAAAGCTTATGTTGTTTAT -p CAGAAGCTTATGTTATTCAT -p CAGAAGCTTATGTTATTTAT -p CAGAAGCTTATGTTGTTCAT -p CAGAAGCTTATGTTGTTTAT -p ATAAACAACATAAGCTTCTGA -p ATAAATAACATAAGCTTCTGA -p ATAAACAACATGAGCTTCTGA -p ATAAATAACATGAGCTTCTGA -p ATAAACAATATAAGCTTCTGA -p ATAAATAATATAAGCTTCTGA -p ATAAACAATATGAGCTTCTGA -p ATAAATAATATGAGCTTCTGA -p ATGAACAACATAAGCTTCTGA -p ATGAATAACATAAGCTTCTGA -p ATGAACAACATGAGCTTCTGA -p ATGAATAACATGAGCTTCTGA -p ATGAACAATATAAGCTTCTGA -p ATGAATAATATAAGCTTCTGA -p ATGAACAATATGAGCTTCTGA -p ATGAATAATATGAGCTTCTGA -m 1 | perl -0076 -ne 'chomp;($h,@S)=split/\n/;$s=join("",@S);print"$h\t$s\n"unless(!$h)' | sed 's/_/ /g' | awk -F " " -v a="$LENGTH" -v b="$COVERAGE" '$4>=a && $6>=b' |  sed 's/ /_/g'  | sort -r -t_ -nk6 | head -${TOP} | sed 's/\t/\n/g' | sed "s/NODE/\>${SAMPLE}_NODE/g" > ../seqkit/coi/${SAMPLE}_coi.fasta
done
chmod -R g+w ../seqkit

#count matches
cat ../seqkit/coi/*coi.fasta | grep ">" | wc -l
#122 # -m 3 with full length-4 (all 4 primers)
#158 # -m 3 with full length-4 and -2 more bp from VertCOI_7194_F
#144 # with full length-4 and -2 more bp from VertCOI_7194_F
#22  # with full length-4 Mod primers only
#13  # -m 3 with full length-4 Mod primers only
#188 # -m 3 with 1/2/2 fewer bp from Mod_F/Mod_R/VertCOI_7194_F
#188 # -m 3 with Mod only, 1/2 fewer bp from Mod_F/Mod_R
#190 # -m 4 with Mod only, 1/2 fewer bp from Mod_F/Mod_R
# 190 # -m 4 with 1/2/2 fewer bp from Mod_F/Mod_R/VertCOI_7194_F #Note: with -m 4 too many are insects! Should we do a directional filtering to figure out which of the matches from -m 4 are vertebrates? using the top coverage is not working, but the question is if there ARE vertebrate sequences present, they are just not top (#2, #3, etc? what distinguishes them?

# 129 # -m 1, 8 no sig match, 2 mosquitoes --> 119 identified
# 159 # -m 2, 2 no sig match, 40 mosquitoes --> 117 identified
# 185 # -m 3, 90 mosquitoes --> 95 identified

```

## For Ace2: 
Primers for ace2: 	ace2-F1457 GAGGAGATGTGGAATCCCAA, 	ace2-B1246s TGGAGCCTCCTCTTCACGG   
Search for primer sequence after trimming 4 bp from start/end, with 3 mismatches allowed  

Add a step to filter to species with internal species diagnostic primers, with 1 mismatch allowed
ACEpip GGAAACAACGACGTATGTACT, trimmed first/last 2, 1 mismatch allowed, alignment shows >2 diagnostic positions
ACEquin CCTTCTTGAATGGCTGTGGCA, trimmed first/last 2, 3 mismatch allowed, alignment shows >4 diagnostic positions

```
bash
module load seqkit
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702/denovo_assembly
mkdir ../seqkit/ace2 ; chmod -R g+w ../seqkit/ace2
mkdir ../seqkit/ACEpip ; chmod -R g+w ../seqkit/ACEpip
mkdir ../seqkit/ACEquin ; chmod -R g+w ../seqkit/ACEquin

LENGTH=150 #change this to set min length
COVERAGE=1 #change this to set min coverage
TOP=1 #change this to set number of top sorted contigs to retain

for SAMPLE in `ls -l | grep -v "total" |  grep -v "fasta" | awk '{print $NF}'`; do
  echo $SAMPLE
  cat ./${SAMPLE}/contigs.fasta | seqkit grep -s -i -p AGATGTGGAATC -p GCCTCCTCTTC -m 3 | perl -0076 -ne 'chomp;($h,@S)=split/\n/;$s=join("",@S);print"$h\t$s\n"unless(!$h)' | sed 's/_/ /g' | awk -F " " -v a="$LENGTH" -v b="$COVERAGE" '$4>=a && $6>=b' |  sed 's/ /_/g'  | sort -r -t_ -nk6 | head -${TOP} | sed 's/\t/\n/g' | sed "s/NODE/\>${SAMPLE}_NODE/g" > ../seqkit/ace2/${SAMPLE}_ace2.fasta
  cat ./${SAMPLE}/contigs.fasta | seqkit grep -s -i -p AGATGTGGAATC -p GCCTCCTCTTC -m 3 | seqkit grep -s -i -p AAACAACGACGTATGTA -m 1 | perl -0076 -ne 'chomp;($h,@S)=split/\n/;$s=join("",@S);print"$h\t$s\n"unless(!$h)' | sed 's/_/ /g' | awk -F " " -v a="$LENGTH" -v b="$COVERAGE" '$4>=a && $6>=b' |  sed 's/ /_/g'  | sort -r -t_ -nk6 | head -${TOP} | sed 's/\t/\n/g' | sed "s/NODE/\>ACEpip_${SAMPLE}_NODE/g" > ../seqkit/ACEpip/${SAMPLE}_ACEpip.fasta
  cat ./${SAMPLE}/contigs.fasta | seqkit grep -s -i -p AGATGTGGAATC -p GCCTCCTCTTC -m 3 | seqkit grep -s -i -p TTCTTGAATGGCTGTGG -m 3 | perl -0076 -ne 'chomp;($h,@S)=split/\n/;$s=join("",@S);print"$h\t$s\n"unless(!$h)' | sed 's/_/ /g' | awk -F " " -v a="$LENGTH" -v b="$COVERAGE" '$4>=a && $6>=b' |  sed 's/ /_/g'  | sort -r -t_ -nk6 | head -${TOP} | sed 's/\t/\n/g' | sed "s/NODE/\>ACEquin_${SAMPLE}_NODE/g" > ../seqkit/ACEquin/${SAMPLE}_ACEquin.fasta
done
chmod -R g+w ../seqkit

#return all ace2 sequences/headers
cat ../seqkit/ace2/*.fasta

#return all ACEpip matching sequences/headers
cat ../seqkit/ACEpip/*fasta
cat ../seqkit/ACEpip/*fasta | grep ">" | awk '{print $1}'

#return all ACEquin matching sequences/headers
cat ../seqkit/ACEquin/*fasta 
cat ../seqkit/ACEquin/*fasta | grep ">" | awk '{print $1}'
```

## For cqm1:  
Primers for cqm1: cqm1-F894 ATGACGGAAGCGTATTCGAG, cqm1-R1834 AAGGTTGATAGCAGCTGCCG    
Search for primer sequence after trimming 4 bp from start/end, with 3 mismatches allowed  

```
bash
module load seqkit
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702/denovo_assembly
mkdir ../seqkit/cqm1 ; chmod -R g+w ../seqkit/cqm1
LENGTH=150 #change this to set min length
COVERAGE=1 #change this to set min coverage
TOP=1 #change this to set number of top sorted contigs to retain

for SAMPLE in `ls -l | grep -v "total" |  grep -v "fasta" | awk '{print $NF}'`; do
  echo $SAMPLE
  cat ./${SAMPLE}/contigs.fasta | seqkit grep -s -i -p CGGAAGCGTATT -p TTGATAGCAGCT -m 3 | perl -0076 -ne 'chomp;($h,@S)=split/\n/;$s=join("",@S);print"$h\t$s\n"unless(!$h)' | sed 's/_/ /g' | awk -F " " -v a="$LENGTH" -v b="$COVERAGE" '$4>=a && $6>=b' |  sed 's/ /_/g'  | sort -r -t_ -nk6 | head -${TOP} | sed 's/\t/\n/g' | sed "s/NODE/\>${SAMPLE}_NODE/g" > ../seqkit/cqm1/${SAMPLE}_cqm1.fasta
done
chmod -R g+w ../seqkit

#return all cqm1 sequences
cat ../seqkit/cqm1/*cqm1.fasta
```


# Step 5: Blast

## Another idea is to use blastn filtered for mammal/birds
This requires running blast on CHPC since there are too many reads to run it online

One last approach we can use to get data from failed COi samples  
a. Blast sorted contigs from spades, length >150 bp, >10X, top 10  
b. Filter results by organism: mammals, birds    
c. Do we get any results that conflict with COi sequences?  

```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/UT-M07101-240702/denovo_assembly
bash
LENGTH=150 #change this to set min length
COVERAGE=10 #change this to set min coverage
NUM=10 #change this to set number of top sorted contigs to retain
for SAMPLE in `ls -l | grep -v "total" |  grep -v "fasta" | awk '{print $NF}'`; do
  echo $SAMPLE
  cat ./${SAMPLE}/contigs.fasta | perl -0076 -ne 'chomp;($h,@S)=split/\n/;$s=join("",@S);print"$h\t$s\n"unless(!$h)' | sed 's/_/ /g' | awk -F " " -v a="$LENGTH" -v b="$COVERAGE" '$4>=a && $6>=b' |  sed 's/ /_/g'  | sort -r -t_ -nk6 | sed 's/\t/\n/g' | sed "s/NODE/\>${SAMPLE}_NODE/g" > ./${SAMPLE}_min${LENGTH}_top${NUM}.fasta
done
chmod -R g+w ../denovo_assembly

#combine all output into one file
cat *top10.fasta > all_min150_top10.fasta    #concatenate top 10 of all sorted matches
# still 18x too many bp to put into one online query

#testing blastn on CHPC
blastn -query B002-UT-M07101-240702_S1_min150_sorted.fasta -db /scratch/general/vast/app-repo/blastdb/DbFiles/v5/nt -out blast_out_B002-UT-M07101-240702_S1_min150_sorted.txt

```

# 09/26/2024 

## Eric and I talked about using a clustering method after filtering the SPAdes contigs for length and coverage  

   a. Clustering all-vs-all: Blast/search all against all, then use e-values as the clustering statistic  
      - FASTA36 as the alignment tool --> Markov clustering in R ("MCL" package)  
      - MMSeqS2 as the alignment and clustering tool https://github.com/soedinglab/MMseqs2  
   b. Clustering all-vs-refs. Same basic options with only reference sequences as the query end of the search. 
   c. Using the primer sequences to find matches.  

## Steps we might need:
  1. FASTQC  
  2. SPAdes de novo assembly  
  3. Filtering contigs with thresholds and add in references for cqm1/ace2  
  4. Clustering and identifying unique clusters  
  5. Assign clusters using the references for cqm1/ace2, but for CO1 we might need to blastn with representatives?  
  6. QC to check for paralogs, haplotype counts, how are SNPs scores/handled, allele balance, this could include mapping reads back to the contigs  

## Assignments:  
Eric will check what we can do with the clustering   
Norah will clean up the readme  
Catch up the other half with our new steps 1 and 2  

