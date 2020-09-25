#! /usr/bin/perl

use strict;
use warnings;

my ($readType) = @ARGV;
if (@ARGV == 0){die "Usage: $0 <PE|SE>\n";}
if(!($readType =~ /PE|SE/)){die "Usage: $0 <PE|SE>\n";}

my $pwd = `pwd`; chomp($pwd);
if(! ($pwd =~ /\/reads/)){die "### ERROR: Please make sure you are located at /path/of/project/reads\n";}

my $checkFastq = `ls *.fastq.gz | wc -l`; chomp($checkFastq);
if($checkFastq == 0){die "### ERROR: No Fastq file in reads/ !\n";}

open OUT, ">input.fofn";

#print OUT "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\tFilePath\n";

if($readType eq "PE"){
	my @fastqs = `ls *_R1_*.fastq.gz | sort | uniq`; chomp(@fastqs);
	foreach my $fastq (@fastqs){
		my @sp = split/_R1/, $fastq;
		my $prefix = $sp[0];
		my $R1_fastq = "$prefix\_R1_001.clean.fastq.gz";
		my $R2_fastq = "$prefix\_R2_001.clean.fastq.gz";
		print OUT "$prefix\t$prefix\treads/$R1_fastq\treads/$R2_fastq\n";
	}
}

if($readType eq "SE"){
	my @fastqs = `ls *.fastq.gz | sort | uniq`; chomp(@fastqs);
	foreach my $fastq (@fastqs){
		my @sp = split/\./, $fastq;
		my $prefix = $sp[0];
		print OUT "$prefix\t$prefix\treads/$fastq\n";
	}
}
close OUT;
