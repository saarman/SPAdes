#!/usr/bin/perl

use strict;
use warnings;

# Output directory
my $output_dir = "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/output";

# Path to MMseqs2 binary
my $mmseqs = "/uufs/chpc.utah.edu/sys/installdir/r8/mmseqs2/oct24/bin/mmseqs";  # Update this path to the actual location of MMseqs2 binary, you can get this by loading module and running: which mmseqs

# Path to the input fasta file
my $fasta = "/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/MMseqs2/input/all.fasta";  # Update this path to your actual input file

# Extract the identifier from the filename
$fasta =~ m/([A-Za-z_\-0-9]+)\.fasta$/ or die "failed match for file $fasta\n";
my $ind = $1;  # Store the identifier in $ind

# Create MMseqs2 database
my $cmd_createdb = "$mmseqs createdb $fasta ${output_dir}/${ind}_DB";
system($cmd_createdb) == 0 or die "system $cmd_createdb failed: $?";

# Run MMseqs2 linclust
my $cmd_linclust = "$mmseqs linclust ${output_dir}/${ind}_DB ${output_dir}/${ind}_DB_lin_clu /scratch/general/vast/u6036559/spades_tmp";
system($cmd_linclust) == 0 or die "system $cmd_linclust failed: $?";

# Create TSV file with cluster information
my $cmd_createtsv = "$mmseqs createtsv ${output_dir}/${ind}_DB ${output_dir}/${ind}_DB ${output_dir}/${ind}_DB_lin_clu ${output_dir}/${ind}_DB_lin_clu.tsv";
system($cmd_createtsv) == 0 or die "system $cmd_createtsv failed: $?";

print "Clustering completed for $ind\n";
