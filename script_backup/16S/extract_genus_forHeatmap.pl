#!/usr/bin/env perl

use warnings;
use strict;

open IN, "<heatmap.txt";
open OUT, ">heatmap.genus.txt";

while(<IN>){
	chomp;
	my @sp = split/\t/, $_;
	my $taxon = $sp[0];
	if ($taxon =~ /g__(.+)$/){
		my $genus = $1;
		print OUT "$genus\t";
	} elsif ($taxon =~ /Taxon/) {
		print OUT "$taxon\t";
	}
	for(my $i=1; $i<=$#sp; $i++){
		print OUT "$sp[$i]\t";
	}
	print OUT "\n";
}
close IN;
close OUT;
