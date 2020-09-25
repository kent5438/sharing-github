#! /usr/bin/env perl

use strict;

my ($txt) = @ARGV;

open IN, "<$txt";
open OUT, ">Mapping_$txt";

my %hash;

print OUT "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\tFilePath\n";
while(<IN>){
	chomp;
	my @sp = split/:/, $_;
	my ($groups, $samples) = @sp;
	my @sample_list = split/,/, $samples;
	foreach my $sample (@sample_list){
		print OUT "$sample\t\t\t$groups\t/export/EC1680U/kentchen/16S/MS18182-MS19042-MS19182.16S/data/$sample.R1.clean.fastq,/export/EC1680U/kentchen/16S/MS18182-MS19042-MS19182.16S/data/$sample.R2.clean.fastq\n";
	my $check_R1 = `ls /export/EC1680U/kentchen/16S/MS18182-MS19042-MS19182.16S/data/$sample.R1.clean.fastq | wc -l`;
	chomp($check_R1);
	if($check_R1 != 0){print "$sample.R1.clean.fastq: *** ok, existed ***\n";}
	else {print "### ERROR: $sample.R1.clean.fastq: *** failed, unmatched ***\n";}
	my $check_R2 = `ls /export/EC1680U/kentchen/16S/MS18182-MS19042-MS19182.16S/data/$sample.R2.clean.fastq | wc -l`;
	chomp($check_R2);
	if($check_R1 != 0){print "$sample.R1.clean.fastq: *** ok, existed ***\n";}
	else {print "### ERROR: $sample.R2.clean.fastq: *** failed, unmatched ***\n";}
	}
}
