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

# Run MMseqs2 cluster with specified number of threads and todo add --min-seq-id .9
my $cmd_cluster = "$mmseqs cluster ${output_dir}/${ind}_DB ${output_dir}/${ind}_DB_clu /scratch/general/vast/u6036559/spades_tmp --threads $threads";
system($cmd_cluster) == 0 or die "system $cmd_cluster failed: $?";

# Create TSV file with cluster information
my $cmd_createtsv = "$mmseqs createtsv ${output_dir}/${ind}_DB ${output_dir}/${ind}_DB ${output_dir}/${ind}_DB_clu ${output_dir}/${ind}_DB_clu.tsv";
system($cmd_createtsv) == 0 or die "system $cmd_createtsv failed: $?";

print "Clustering completed for $ind\n";

$pm->finish;  # End the child process

$pm->wait_all_children;  # Wait for all child processes to finish
