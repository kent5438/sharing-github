#! /usr/bin/env perl

use strict;
use warnings;

my ($path) = @ARGV;

#open IN, "<$path/meme.xml";
#open OUT, ">$path/tmp1.txt";
#while(<IN>){
#	chomp;
#	if(/<sequence id=(.+)name=\"(.+)\" length=\"([0-9]+)\"/){
#		print OUT "$3\n"; 
#	}
#}
#close IN;
#close OUT;

open IN, "<$path/meme.txt";
open OUT, ">$path/motif_summary.txt";
while(<IN>){
	chomp;
	if(/(M02662:177:000000000-CPB)\s+([0-9]+)\s+((.+)e-[0-9]+)\s+(.+)/){
		print OUT "$1\t$3\t$5\t";
		my $CA_rep_count = `echo $5 | sed 's/CA/z/g' | grep -oP 'z{2,}' | grep -o 'z' | wc -l`; chomp $CA_rep_count;
		my $TG_rep_count = `echo $5 | sed 's/TG/x/g' | grep -oP 'x{2,}' | grep -o 'x' | wc -l`; chomp $TG_rep_count;
		my $rep_count = $CA_rep_count + $TG_rep_count;
		print OUT "$rep_count\t";
		my $count = `echo $5 | sed 's/ //g' | grep -o '.' | wc -l`; chomp $count;
		print OUT "$count\n";
	}
	if(/Motif (.+) MEME-1 regular expression/){
		print OUT "Motif Consensus Sequence\t\t$1\n";
	}
}
close IN;
close OUT;



system("sed -i '1i seq_id\tpvalue\tmotif_seqs\tconsecutive_CA-TG_count\tlength' $path/motif_summary.txt");
