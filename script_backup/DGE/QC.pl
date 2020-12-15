#! /usr/bin/perl

use warnings;
use strict;

my $USER = `id -u -n`;
chomp($USER);

my @reads = `ls reads/*.fastq.gz`;
chomp(@reads);

if(scalar @reads == 0){die "### ERROR: no read existed in 'reads/' folder";}

if(! -e "fastqc_result"){mkdir("fastqc_result");}
my @fastqc_content = `ls fastqc_result`;
chomp(@fastqc_content);
if(scalar @fastqc_content == 0){
	my $fastqc = "/export/EC1680U/software/FastQC/fastqc";
	system("$fastqc -t 32 -o fastqc_result @reads");
}

my $multiqc = "/export/EC1680U/software/anaconda2/bin/multiqc";
if(! -e "multiqc_report.html"){
	my $multiqc_cmd = "$multiqc -f ".
	"-c /export/EC1680U/Pipeline/DGE/multiqc_config.yaml ".
	"fastqc_result output/1.CalculateExpression";
	system("$multiqc_cmd");
}

my $html_path = "output/5.report/html";
if(! -e $html_path){die "### ERROR: make sure $html_path exist!\n";}
system("rsync -avq multiqc_* $html_path");
