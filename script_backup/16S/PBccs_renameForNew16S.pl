#! /usr/bin/env perl

use warnings;
use strict;

my $pwd = `pwd`; chomp($pwd);

my ($barcode) = @ARGV;
if(! $barcode){die "\n### Usage: $0 <PS17015-16S-barcode_1.txt>\n\n";}

my $check_ccs = `ls *.fastq | wc -l`; chomp($check_ccs);
if($check_ccs != 0){
	my @ccses = `ls *.fastq | sort -V`; chomp(@ccses);
	my @samples = `cat $barcode | grep -v sampleName | cut -f4`; chomp(@samples);

	my $i=0;
	foreach my $ccs (@ccses){
		system("mv $ccs $samples[$i].fastq");
		$i++;
	}
} else {die "### Please make sure CCS fastq are located in this folder or the filenames are corrected! \n";}

