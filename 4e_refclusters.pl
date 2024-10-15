#!/usr/bin/perl

use strict;
use warnings;

# Number of threads to use with MMseqs2
my $threads = 20;  # Adjust this based on your SLURM settings and system

# Paths to the input, output, and MMseqs2 binary
my $output_dir = "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output";
my $mmseqs = "/uufs/chpc.utah.edu/sys/installdir/r8/mmseqs2/oct24/bin/mmseqs";

# Step 1: Generate a pseudo-FASTA file for the clusters matching refs
my $resultsDB = "$output_dir/searchRepsDB"; # Search results DB
my $queryDB = "$output_dir/DB_refs"; # DB for query
my $targetDB = "$output_dir/DB_rep"; # DB for subject, already exists!!!
my $db_seq= "$output_dir/DB_fsa"; # Output for the database
my $out_fasta= "$output_dir/ref_clu.fasta"; # Output for the fasta
my $createseqfiledb_cmd = "$mmseqs createseqfiledb $resultsDB $resultsDB $targetDB"; # Create database
my $result2flat_cmd = "$mmseqs result2flat $resultsDB $resultsDB $targetDB $out_fasta"; # Create fasta
