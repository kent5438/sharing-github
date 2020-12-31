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

open OUT, ">../input.fofn";

#print OUT "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\tFilePath\n";

if($readType eq "PE"){
	my @fastqs = `ls *_R1*.fastq.gz | awk -F'_R1' '{print \$1}'`; chomp(@fastqs);
	foreach my $fastq (@fastqs){
		my @sp = split/_/, $fastq;
		my $prefix = $sp[0];
		my $R1_fastq = "$prefix\_R1.clean.fastq.gz";
#		system("mv $fastq\_R1.clean.fastq.gz $R1_fastq");
		my $R2_fastq = "$prefix\_R2.clean.fastq.gz";
#		system("mv $fastq\_R2.clean.fastq.gz $R2_fastq");
		print OUT "$prefix\t$prefix\treads/$R1_fastq\treads/$R2_fastq\n";
	}
}

if($readType eq "SE"){
	my @fastqs = `ls *.fastq.gz | awk -F'_R1' '{print \$1}'`; chomp(@fastqs);
	foreach my $fastq (@fastqs){
		my @sp = split/\_/, $fastq;
		my $prefix = $sp[0];
		my $R1_fastq = "$prefix\_R1.clean.fastq.gz";
		system("mv $fastq\_R1.clean.fastq.gz $R1_fastq");
		print OUT "$prefix\t$prefix\treads/$R1_fastq\n";
	}
}
close OUT;
