#! /usr/bin/perl

use strict;

open IN, "<all_species.txt";
open OUT, ">species_hash.txt";

my %hash;
while(<IN>){
	chomp;
	$hash{$_}++;
}
foreach my $species (sort keys %hash){
	print OUT "$species\t$hash{$species}\n";
}
close IN;
close OUT;
