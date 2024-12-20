#!/usr/bin/perl

use strict;
use warnings;

# Number of threads to use with MMseqs2
my $threads = 20;  # Adjust this based on your SLURM settings and system

# Paths to the output directory and MMseqs2 binary and temp dir
my $output_dir = "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output";
my $mmseqs = "/uufs/chpc.utah.edu/sys/installdir/r8/mmseqs2/oct24/bin/mmseqs";
my $tmp_dir = "/scratch/general/vast/u6036559/spades_tmp/";    # Temporary directory for MMseqs2

# Paths to inputs query and target
my $ref_fasta = "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input/mmRefs.fasta"; # Query
my $rep_fasta = "$output_dir/clu_rep.fasta"; # Target

# Step 1: Generate a m8 from easy-search with ref seqs against representatives for each cluster
my $out_m8 = "$output_dir/easyReps.m8"; # Search m8 output
my $easy_search_cmd = "$mmseqs easy-search $ref_fasta $rep_fasta $out_m8 $tmp_dir --search-type 3"; # Search

# Step 2: Execute the commands
print "Running MMseqs2 easy-search...\n";
system($easy_search_cmd) == 0 or die "MMseqs2 easy-search command failed: $?";

print "MMseqs2 search tasks completed successfully.\n";
