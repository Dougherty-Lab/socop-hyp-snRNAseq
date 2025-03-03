#!/usr/bin/perl

use strict;
use warnings;

my $usage = "\n\n\tcreateI1.pl FASTQ-file\n\n";

die $usage unless ( @ARGV==1);

my $file = $ARGV[0];

open (IN, "gunzip -c $file | ") or die "Can't open gzipped file: $file\n";

while (<IN>){

    if ( /^\@.*:[ACGTN]{6,20}\+([ACGTN]{6,20})/ ){

	my $index = $1;
	my $len   = length($index);
	my $qual  = "?" x $len;
	print $_, $index, "\n+\n", $qual, "\n";
    }
}

# Note: output I1 is not gzipped
# After running this script, run: 
# gzip <*_I1.fastq>