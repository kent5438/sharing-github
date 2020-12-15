#! /usr/bin/env perl
# Author: Kent Chen
# Date: 2020/07/02
# Application: Eukaryotic genome de-novo annotation pipeline (without RNA)
# version: 0.1

use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use Time::Seconds;
use lib '/export/EC1680U/kentchen/lib';
use Time qw(CurrTime);


my $pwd = `pwd -L`; chomp($pwd);
my $USER = `id -u -n`; chomp($USER);
my $tab2xls = "/export/EC1680U/perl/bin/tab2xls.pl";
my $tab2xlsx = "/export/EC1680U/perl/bin/tab2xlsx.pl";
my ($input, $output, $species);

my ($help, $prefix, $code);
GetOptions(
	"help!"	=> \$help,
	"input=s"	=> \$input,
	"output=s"	=> \$output,
	"species=s"	=> \$species,
);

if($help){
	die "\n### USAGE ###
* Used for Eukaryotic genome de-novo annotation pipeline (without RNA)
- Funannotate (main package)
- redmask.py (repeat mask tool)

E.g.)
% $0 -i <scaffolds.polished.fasta> -o funannotate -s ustilago

--input/-i : final assembled genome fasta
--output/-o : output path
--species/-s : used for training set model organism
"
} elsif (! $input){
	die "\n### ERROR: The input genome fasta is necessary, e.g. scaffolds.polished.fasta\n";
} elsif (! $output){
	die "\n### ERROR: output path is necessary e.g. funannotate\n";
} elsif (! $species){
	die "\n### ERROR: training species is necessary, e.g. ustilago\n
try: % funannotate species";
} else {
	open OUT, ">>CMD_history.txt";
	print OUT CurrTime() . "$0 -i $input -o $output -s $species\n";
	close OUT;
}



if(! -e "funannotate"){
	mkdir("funannotate");
}
chdir("funannotate");




print CurrTime() . "Complete!\n";
