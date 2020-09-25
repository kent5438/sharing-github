#! /usr/bin/env perl

use strict;
use warnings;

my %hash;
my (@keys, @values);
my $pwd = `pwd -L`; chomp($pwd);
my $r1_fq = `ls *R1* | wc -l`; chomp($r1_fq);
my $r2_fq = `ls *R2* | wc -l`; chomp($r2_fq);
my $seqtk = '/export/EC1680U/kentchen/miniconda3/bin/seqtk';

my @fq_counts = `grep -c '^\@M02662' *?R?.clean.fastq*`; chomp(@fq_counts);
foreach my $fq_count (@fq_counts){
	my @sp = split/:/, $fq_count;
	my ($fq, $count) = @sp;
	$hash{$fq} = $count;
	if($count >= 100000 && $count < 120000){
		my $cmd = "$seqtk sample -s100 $fq 0.7 > $fq.sub;";
		$cmd .= "mv $fq $fq\.raw;";
		$cmd .= "mv $fq\.sub $fq";
		system("$cmd");
	} 
	if($count >= 120000 && $count < 140000){
		my $cmd = "$seqtk sample -s100 $fq 0.6 > $fq.sub;";
		$cmd .= "mv $fq $fq\.raw;";
		$cmd .= "mv $fq\.sub $fq";
		system("$cmd");
	} 
	if($count >=140000 && $count < 160000){
		my $cmd = "$seqtk sample -s100 $fq 0.5 > $fq.sub;";
		$cmd .= "mv $fq $fq\.raw;";
		$cmd .= "mv $fq\.sub $fq";
		system("$cmd");
	}
	if($count >= 160000){
		my $cmd = "$seqtk sample -s100 $fq 0.5 > $fq.sub;";
		$cmd .= "mv $fq $fq\.raw;";
		$cmd .= "mv $fq\.sub $fq";
		system("$cmd");
	}
}

