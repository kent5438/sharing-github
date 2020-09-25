#! /usr/bin/env perl

use strict;
use warnings;

my %hash;
my (@keys, @values);
my $pwd = `pwd -L`; chomp($pwd);
my $seqtk = '/export/EC1680U/kentchen/miniconda3/bin/seqtk';

my @fq_counts = `echo \"\$(ls *.trimmed.fa):\$(grep -c '>' *.trimmed.fa)\"`; chomp(@fq_counts);
foreach my $fq_count (@fq_counts){
	my @sp = split/:/, $fq_count;
	my ($fq, $count) = @sp;
	$hash{$fq} = $count;
	if($count >= 200000 && $count < 250000){
		my $cmd = "$seqtk sample -s100 $fq 0.75 > $fq.sub;";
		$cmd .= "mv $fq $fq\.raw;";
		$cmd .= "mv $fq\.sub $fq";
		system("$cmd");
	} 
	if($count >= 250000 && $count < 300000){
		my $cmd = "$seqtk sample -s100 $fq 0.6 > $fq.sub;";
		$cmd .= "mv $fq $fq\.raw;";
		$cmd .= "mv $fq\.sub $fq";
		system("$cmd");
	} 
	if($count >= 300000){
		my $cmd = "$seqtk sample -s100 $fq 0.5 > $fq.sub;";
		$cmd .= "mv $fq $fq\.raw;";
		$cmd .= "mv $fq\.sub $fq";
		system("$cmd");
	}
}

