#!/usr/bin/perl

use strict;
use warnings;

my %hash;
my ($file) = @ARGV;
open IN, "<$file";
while(<IN>){
	chomp;
	my @sp = split/\t/, $_;
	my ($contig, $cov) = ($sp[0], $sp[3]);
	$hash{$contig}{cov} += $cov;
	$hash{$contig}{count}++;
}

foreach my $contig (sort keys %hash){
	my $cov = $hash{$contig}{cov};
	my $count = $hash{$contig}{count};
	my $avg = $cov / $count;
	print "$contig\t";
	printf "%.2f\n", $avg;
}
close IN;
