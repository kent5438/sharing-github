#! /usr/bin/env perl
### usage: perl bamCount.pl 5699.bam-readCount.txt > 5699.baseCount.txt

use strict;
use warnings;

my ($file) = @ARGV;
open IN, "<$file";
while(<IN>){
	chomp;
	my @sp = split/\t/, $_;
	print "$sp[0]\t$sp[1]\t$sp[2]\t$sp[3]\t";
	my @sp5 = split/:/, $sp[5];
	my @sp6 = split/:/, $sp[6];
	my @sp7 = split/:/, $sp[7];
	my @sp8 = split/:/, $sp[8];
	print "$sp5[0]: $sp5[1]\t";
	print "$sp6[0]: $sp6[1]\t";
	print "$sp7[0]: $sp7[1]\t";
	print "$sp8[0]: $sp8[1]\t";
	
	if($sp[3] == 0){
		print "$sp5[0]: 0%\t";
		print "$sp6[0]: 0%\t";
		print "$sp7[0]: 0%\t";
		print "$sp8[0]: 0%\n";
	} else {
		my $sp5_pct = $sp5[1] / $sp[3] * 100;
		my $sp6_pct = $sp6[1] / $sp[3] * 100;
		my $sp7_pct = $sp7[1] / $sp[3] * 100;
		my $sp8_pct = $sp8[1] / $sp[3] * 100;
		print "$sp5[0]: $sp5_pct%\t";
		print "$sp6[0]: $sp6_pct%\t";
		print "$sp7[0]: $sp7_pct%\t";
		print "$sp8[0]: $sp8_pct%\n";
	}
}
close IN;
