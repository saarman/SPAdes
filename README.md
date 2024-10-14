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

# Step 4: Clustering with MMseqs2

Example of the raw code:
```
# Load modules
module load mmseqs2/oct24  # change to module name

# Check if working by loading help menu
mmseqs --help

## convert fasta to DB
mmseqs createdb ./input/file.fasta DB

## cluster with cluster or linclust (linclust run time scales linearly but is slightly less accurate)
# mmseqs cluster DB DB_clu tmp
# clusters can be set by adding --min-seq-id .9
# make sure tmp folder exists
mmseqs linclust DB DB_lin_clu /scratch/general/vast/u6036559/spades_tmp

## outputting files
# TSV file with representative cluster sequences on left and all members of the cluster on right
mmseqs createtsv DB DB DB_lin_clu DB_lin_clu.tsv
```

## Example of the sbatch script: 4a_MMseqs2.slurm
```
#!/bin/sh
#SBATCH --time=336:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20          # same as $max set in ForkManager
#SBATCH --account=saarman-np
#SBATCH --partition=saarman-shared-np   
#SBATCH --job-name=MMseqs2_try1
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=norah.saarman@usu.edu

# Load modules
module load mmseqs2/oct24  # change to module name

# Change to the directory where the input data is located
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input

# Concatenate all of the input into one file
rm all.fasta; cat *.fasta > all.fasta
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/

# Run the Perl script with the input files
perl /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/scripts/4a_MMseqs2.pl all.fasta

# Permissions
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/
```

## Example of the perl script: 4a_MMseqs2.pl
```
#!/usr/bin/perl

use strict;
use warnings;
use Parallel::ForkManager;

my $max = 1;  # Set the maximum number of parallel processes to 1 since we are managing threads within MMseqs2
my $pm = Parallel::ForkManager->new($max);  # Create a new Parallel::ForkManager object with the specified maximum

# Output directory
my $output_dir = "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output";

# Path to MMseqs2 binary
my $mmseqs = "/uufs/chpc.utah.edu/sys/installdir/r8/mmseqs2/oct24/bin/mmseqs";  # Update this path to the actual location of MMseqs2 binary

# Path to the input fasta file
my $fasta = "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input/all.fasta";  # Update this path to your actual input file

# Number of threads to use for MMseqs2
my $threads = 20;  # Adjust this number based on your available resources

$pm->start and next if $pm->start;  # Fork a new process

# Extract the identifier from the filename
$fasta =~ m/([A-Za-z_\-0-9]+)\.fasta$/ or die "failed match for file $fasta\n";
my $ind = $1;  # Store the identifier in $ind

# Create MMseqs2 database
my $cmd_createdb = "$mmseqs createdb $fasta ${output_dir}/${ind}_DB";
system($cmd_createdb) == 0 or die "system $cmd_createdb failed: $?";

# Run MMseqs2 linclust with specified number of threads
my $cmd_linclust = "$mmseqs linclust ${output_dir}/${ind}_DB ${output_dir}/${ind}_DB_lin_clu /scratch/general/vast/u6036559/spades_tmp --threads $threads";
system($cmd_linclust) == 0 or die "system $cmd_linclust failed: $?";

# Create TSV file with cluster information
my $cmd_createtsv = "$mmseqs createtsv ${output_dir}/${ind}_DB ${output_dir}/${ind}_DB ${output_dir}/${ind}_DB_lin_clu ${output_dir}/${ind}_DB_lin_clu.tsv";
system($cmd_createtsv) == 0 or die "system $cmd_createtsv failed: $?";

print "Clustering completed for $ind\n";

$pm->finish;  # End the child process

$pm->wait_all_children;  # Wait for all child processes to finish

```

Before running, I need to make these files on github, and then use git to clone
```
# Just once:
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2
git clone https://github.com/saarman/SPAdes scripts
```

## Every time need to Pull and use sbatch to run:
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/scripts
git pull
sbatch 4a_MMseqs2.slurm
```
Output on tsv format: /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output/all_DB_lin_clu.tsv

## Visualize the results R script named 4b_plot_clusters.R
```
# Load required libraries
library(ggplot2)
library(data.table)
library(igraph)
library(ggraph)

# Define the input file path
input_file <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output/all_DB_lin_clu.tsv"

# Define the output directory
output_dir <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output/"  # Update this path to your desired output directory

# Read the TSV file
data <- fread(input_file, header = FALSE)

# Assign column names (assuming the format is representative sequence and cluster members)
colnames(data) <- c("Representative", "ClusterMember")

# Create an edge list for the network graph
edge_list <- data[, .(Representative, ClusterMember)]

# Create an igraph object from the edge list
network_graph <- graph_from_data_frame(edge_list, directed = FALSE)

# Generate layout coordinates
layout <- create_layout(network_graph, layout = "fr")

# Identify representative nodes
representative_nodes <- unique(edge_list$Representative)

# Plot the network with only representative nodes labeled
p <- ggraph(layout) +  # Use the layout with coordinates
  geom_edge_link(width = 1) +  # Edges without color
  geom_node_point(size = 2) +  # Size of nodes
  geom_node_text(aes(label = ifelse(name %in% representative_nodes, name, "")), vjust = 1, hjust = 1) +  # Label only representative nodes
  theme_minimal() +  # Minimal theme
  ggtitle("Clustering Network")

# Save the plot to a file in the specified output directory
output_file <- file.path(output_dir, "network_plot.png")
ggsave(output_file, plot = p, width = 10, height = 7)
```
... to save network_plot.png in /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output/
  
Pull and run R script
```
cd /uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/scripts
git pull
module load R
Rscript 4b_plot_clusters.R
```
