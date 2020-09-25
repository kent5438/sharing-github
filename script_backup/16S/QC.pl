#! /usr/bin/env perl

use warnings;
use strict;

my $USER = `id -u -n`;
chomp($USER);

my @reads = `ls 0_cleanFastq/*.fastq`; chomp(@reads);
my $read_R1R2_check = `ls 0_cleanFastq/*.fastq | grep -c 'R1'`; chomp($read_R1R2_check);

if(scalar @reads == 0){die "### ERROR: no read existed in 'reads/' folder";}


if(! -e "fastqc_result"){mkdir("fastqc_result");}
my @fastqc_content = `ls fastqc_result`;
chomp(@fastqc_content);
if(scalar @fastqc_content == 0){
	my $fastqc = "/export/EC1680U/software/anaconda2/bin/fastqc";
	system("$fastqc -t 32 -o fastqc_result @reads");
}

my $multiqc = "/export/EC1680U/software/anaconda2/bin/multiqc";
if(! -e "multiqc_report.html"){
	if ($read_R1R2_check != 0) {
		my $multiqc_cmd = "$multiqc ".
	    "-c /export/EC1680U/Pipeline/16S/MiSeq-16S_multiqc_config.yaml ".
	    "fastqc_result fastqc_result";
	    system("$multiqc_cmd");
	} else {
		my $multiqc_cmd = "$multiqc ".
	    "-c /export/EC1680U/Pipeline/16S/PB-16S_multiqc_config.yaml ".
	    "fastqc_result fastqc_result";
	    system("$multiqc_cmd");
	}
}

my $html_path = "report/html";
if(! -e $html_path){die "### ERROR: make sure $html_path exist!\n";}
system("rsync -avq multiqc_* $html_path");

