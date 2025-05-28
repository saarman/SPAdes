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

## Using Git in CHPC (do this once)
1. cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles
2. git clone https://github.com/saarman/SPAdes SPAdes_scripts

## Need to run everytime 
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts
git pull

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

# Step 4: Search/BLAST with MMseqs2
### Step 4a: MMseqs2 Search with BOLD as reference to Get Top Match per ref
4a_MMseqs2_easyBOLD.slurm
```
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
```
## Step 4b: Filter results
Still need to update Step 4b: Filter for e-val threshold, percent identity, alignment length

These are quality control steps to make sure the hit is a good COI hit, not just the “least bad” one:

E-value threshold (e.g., E < 1e-5 or stricter like 1e-20). Helps remove weak, possibly spurious matches. Even the “top hit” could be a poor match if the contig is junk
Percent identity (e.g., ≥ 85–90%) Indicates how similar your contig is to a known COI. Useful to exclude poor alignments, pseudogenes, or sequencing errors
Alignment length or coverage. Short matches may not be biologically meaningful. Filtering for alignment length >100 bp (or % of COI length) can help


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
sbatch 4a_MMseqs2_easyBOLD.slurm
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







