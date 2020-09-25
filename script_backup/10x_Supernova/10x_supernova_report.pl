#! /usr/bin/env perl

### Build report by R markdown
### Author: KentChen
### Date: 2017-12-28

use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use Time::Seconds;

### set my defaults
my $help;
my $demulx = "mkfastq";
my $asm = "assembly";
my $pwd = `pwd -L`; chomp($pwd);

GetOptions(
	"help!"	=> \$help,
	"demultiplex=s"	=> \$demulx,
	"assembly=s"	=> \$asm,
);

if($help){
	die "\n### USAGE ###
* Used for 10x Supernova Report

E.g.)
% $0 -d mkfastq -a assembly

--demultiplex/-d : demultiplexing folder generating from 'mkfastq' (or read folder)
--assembly/-a : assembly folder generating from 'run'
"
} else {
	print "\n### Input Options ###\n";
	print "-d): $demulx\n";
	print "-a): $asm\n";
	print "% perl $0 -d $demulx -a $asm\n\n";
	open OUT, ">>CMD_history.txt";
	print OUT CurrTime() . "perl $0 -d $demulx -a $asm\n # suggested absolutely path\n";
	close OUT;
}


### Copy report_template to current dir)
print CurrTime() . "Copy report_template to current dir...\n";
my $report_template = "/export/EC1680U/Pipeline/10x_Supernova/report_template";
if(! -e "$report_template"){die "\n### ERROR: Please make sure '$report_template' existed!\n\n";}
system("rsync -avq $report_template .");


### Transfer demultiplex result to report
print CurrTime() . "Transfer demultiplex result to report...\n";
if(! -e $demulx){die "\n### ERROR: Please make sure $demulx/ existed!\n\n";}
if(-e "$demulx/Reports/html"){
	system("rsync -avq $demulx/Reports/html report_template/files/Demultiplex_report");
} else {
	print "\n### WARNING: 'Your fastq might be outsource's data\n";
	print "--- do fastqc/multiqc jobs ---\n";
	if(! -e "fastqc_result"){
		mkdir("fastqc_result");
		system("fastqc -t 8 -o fastqc_result $demulx/*.fastq.gz");
	}
	system("multiqc --force -c /export/EC1680U/Pipeline/10x_Supernova/multiqc_fastqc.yaml -f fastqc_result/ -o report_template/files");
}

### Transfer assembly result to report
print CurrTime() . "Transfer assembly result to report...\n";
if(! -e $asm){die "\n### ERROR: Please make sure $asm/ existed!\n\n";}

if(! -e "$asm/outs/report.txt"){die "\n### ERROR: 'supernova run' completed smoothly?\n\n";}
system("rsync -avq $asm/outs/report.txt report_template/files/Assembly_report.txt");

if(! -d "$asm/quast_results"){die "\n### ERROR: Please make sure 'QUAST' has been comleted!!!\n\n";}
system("rsync -avq $asm/quast_results report_template/files");

### Build R markdown Script
print CurrTime() . "Build Rmarkdown-Style Report...\n";
my $rscript = "/export/EC1680U/software/anaconda2/bin/Rscript";
system("$rscript /export/EC1680U/Pipeline/10x_Supernova/build_site.R");


### Remove report_template & Transfer 4 types genome into report/files
print CurrTime() . "Rsync report to root folder & remove report_template...\n";
system("rsync -avq report_template/report .");
system("rm -rf report_template");

print CurrTime() . "Transfer 4 types genome into report/files...\n";
my $num_of_asm = `ls $asm/*.fasta.gz | wc -l`; chomp($num_of_asm);
if($num_of_asm != 5){die "\n### ERROR: 'supernova mkoutput' completed smoothly?\n\n";}
#system("ln -s $pwd/$asm/*fasta.gz $pwd/report/files");

print CurrTime() . "Complete!\n";

### Subroutine
sub CurrTime {
    my $t = localtime;
    my $time = "[".$t->hms."]*** ";
    return "$time";
}
