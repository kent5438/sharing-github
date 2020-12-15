#! /usr/bin/perl

use warnings;
use strict;

my $pwd = `pwd`; chomp($pwd);
my $data_prefix = `basename $pwd`; chomp($data_prefix);

if($data_prefix ne 'data'){
	die "### ERROR: Please make sure you are located in */data";
}

my ($barcode) = @ARGV;
if(! $barcode){die "\n### Usage: 
* cd data/
* $0 PS17015-16S-barcode_1.txt\n\n";}

my $check_ccs = `ls Ibc*.fastq | wc -l`; chomp($check_ccs);
if($check_ccs != 0){
	my @ccses = `ls Ibc*.fastq | sort -V`; chomp(@ccses);
	my @samples = `cat $barcode | grep -v sampleName | cut -f4`; chomp(@samples);

	my $i=0;
	foreach my $ccs (@ccses){
		system("mv $ccs $samples[$i].clean.fastq");
		$i++;
	}
}

system("cut -f4 $barcode | grep -v sampleName >> ../SampleName_16S.txt");

if((! -e "../1_joinPairs") or (! -e "../0_cleanFastq")){
	system("mkdir ../0_cleanFastq ../1_joinPairs");
}
system("ln -sf $pwd/*.clean.fastq $pwd/../0_cleanFastq");


my @prefixes = `ls *.clean.fastq | cut -d'.' -f1`; chomp(@prefixes);
foreach my $prefix (@prefixes){
	system("ln -sf $pwd/$prefix.clean.fastq $pwd/../1_joinPairs/$prefix.assembled.fastq");
}
