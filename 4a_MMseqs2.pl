#!/usr/bin/perl

use strict;
use warnings;

# Number of threads to use with MMseqs2
my $threads = 20;  # Adjust this based on your SLURM settings and system

# Paths to the input, output, and MMseqs2 binary
my $output_dir = "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output";
my $mmseqs = "/uufs/chpc.utah.edu/sys/installdir/r8/mmseqs2/oct24/bin/mmseqs";
my $fasta = "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input/all.fasta";

# Step 1: Create MMseqs2 database from the input FASTA file
my $db_dir = "$output_dir/DB";  # Directory to store the database
my $createdb_cmd = "$mmseqs createdb $fasta $db_dir";

# Step 2: Cluster sequences with MMseqs2
my $db_clu = "$output_dir/DB_clu";  # Output for the cluster
my $tmp_dir = "/scratch/general/vast/u6036559/spades_tmp/";    # Temporary directory for MMseqs2
my $cluster_cmd = "$mmseqs cluster $db_dir $db_clu $tmp_dir --min-seq-id 0.90 --threads $threads";

# Step 3: Generate a TSV file for the clusters
my $tsv_file = "$output_dir/DB_clu.tsv";
my $createtsv_cmd = "$mmseqs createtsv $db_dir $db_dir $db_clu $tsv_file";

# Step 4: Generate a FASTA file for the clusters
my $db_seq= "$output_dir/DB_fsa";
my $out_fasta= "$output_dir/clu_fasta.fasta";
my $createtsv_cmd = "$mmseqs createseqfiledb $db_dir $db_clu $db_seq; $mmseqs result2flat $db_dir $db_dir $out_fasta";

mmseqs createseqfiledb DB clu clu_seq 
mmseqs result2flat DB DB clu_seq clu_seq.fasta

# Step 5: Execute the commands
print "Running MMseqs2 createdb...\n";
system($createdb_cmd) == 0 or die "MMseqs2 createdb command failed: $?";

print "Running MMseqs2 cluster...\n";
system($cluster_cmd) == 0 or die "MMseqs2 cluster command failed: $?";

print "Running MMseqs2 createtsv...\n";
system($createtsv_cmd) == 0 or die "MMseqs2 createtsv command failed: $?";

print "Running MMseqs2 createtsv...\n";
system($createfasta_cmd) == 0 or die "MMseqs2 createfasta command failed: $?";

print "MMseqs2 tasks completed successfully.\n";
