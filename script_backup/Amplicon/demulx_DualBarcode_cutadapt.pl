#! /usr/bin/env perl

### [Notice]: Input file should be joined fastq file created by PEAR/FLASH

use strict;
use warnings;
use Getopt::Long;

my $pwd = `pwd -L`; chomp($pwd);
my $cutadapt = "/export/EC1680U/software/anaconda2/bin/cutadapt";
my $output = "demultiplex_cutadapt";

if(! -e $output){
	mkdir $output;
}
chdir("$output");

my ($help, $barcode, $fastq);
GetOptions(
	"help!"	=>	\$help,
	"barcode=s"	=>	\$barcode,	# should be 4 columns with header (#dualbarcodes	barcode1	barcode2	sampleName)
	"fastq=s"	=>	\$fastq,	# single-end joined fastq
);

if($help){
	die "\n### USAGE: 
% $0 -b <barcode file> -f <fastq file>

--barcode/-b : should be 4 columns with header (#dualbarcodes    barcode1    barcode2    sampleName)
--fastq/-f : single-end joined fastq
\n";
} elsif(! $barcode){
	die "\n### ERROR: barcode file is required!\n\n";
} elsif(! $fastq){
	die "\n### ERROR: joined fastq file is required";
} else {
	open OUT, ">>CMD_history.txt";
	print OUT "$0 -b $barcode -f $fastq\n";
	close OUT;
}


open IN, "<$pwd/$barcode";
while(<IN>){
	chomp;
	next if $. == 1;
	my @sp = split/\t/, $_;
	my ($demultx, $index1, $index2, $prefix) = @sp;
	my $index1_rev = `echo $index1 | grep '^[ATCG]' | rev | tr ATCG TAGC`; chomp($index1_rev);
	my $index2_rev = `echo $index2 | grep '^[ATCG]' | rev | tr ATCG TAGC`; chomp($index2_rev);
	system("$cutadapt -a $index1...$index2_rev\$ -a $index2...$index1_rev\$ ".
		"--discard-untrimmed --no-trim ".
		"-o $prefix\.clean.fastq ".
		"$pwd/$fastq 2>&1 | tee $prefix.cutadapt.log");
}
close IN;


