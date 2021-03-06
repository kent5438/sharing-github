#! /usr/bin/env perl
# Author: kentchen
# Create date: 2017-06-12
# Last update: 2017-10-02

use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use Time::Seconds;
use lib '/export/EC1680U/kentchen/lib';
use Time qw(CurrTime);


my $USER = `id -u -n`;
chomp($USER);
my $help;
my $input;
GetOptions(
	"help!"	=> \$help,
	"input=s"	=> \$input,
);
#my $fastq_path = "/export/EC2480U/QCResults/fastq";
my $fastq_path = "/mnt/NFS/EC2480U-P/scratch/QCResults/fastq";

if($help){
	die "\n### USAGE ###
* De novo transcriptome assembly by Trinity and following with DGE pipeline

E.g.)
% Trinity_pipe.pl -i input.fofn

-i|--input:
Sample information & Clean read repository path
e.g (gzip also could be used)
	ctrl_1  ctrl_1  reads/ctrl_1.R1.fastq   reads/ctrl_1.R2.fastq
	treat_1 treat_1 reads/treat_1.R1.fastq  reads/treat_1.R2.fastq
	treat_2 treat_2 reads/treat_2.R1.fastq  reads/treat_2.R2.fastq\n\n";
} elsif (! $input){
    die "\n### Error: input file is neccessary\n";
} else {
	open OUT, ">>CMD_history.txt";
	print OUT CurrTime() . "$0 -i $input\n";
	close OUT;
}

### Check reads are matched with input.fofn
# Generate Sample_list.txt
my (@r1s, @r2s);
open my $sample_in, "<$input" or die "### ERROR: please check the '$input' existed or not\n";
my $log = "log";
if(! -e $log){mkdir("$log");}

while(<$sample_in>){
	chomp($_);
	my @sp = split/\t/, $_;
	my ($cond, $sample, $r1, $r2) = @sp;
	
	if(defined ($r1 or $r2)) {
        print "$r1\n$r2\n";
        push(@r1s, $r1);
        push(@r2s, $r2);
    }
	if(undef $r1){
		die "### ERROR: Read 1 is necessary!\n";
	} elsif (undef $r2) {
		print "### WARNINGS: $sample is no Read 2 file, continued with single-end pipeline!\n";
		push(@r1s, $r1);
	}
}
close $sample_in;

### Read quality evaluation by MultiQC
print CurrTime() . "Read Quality Evaluation by MultiQC...\n";
if(! -e "fastqc_result"){mkdir("fastqc_result");}
my $fastqc_content = `ls fastqc_result | wc -l`;
chomp($fastqc_content);
if($fastqc_content == 0){
	my $fastqc_cmd = "fastqc -t 32 -o fastqc_result @r1s @r2s";
    system("$fastqc_cmd");
}

if(! -e "multiqc_fastqc.html"){
	system("/export/EC1680U/perl/bin/Trinity/fastq-stats-k.pipe.sh 150");
}

#my $multiqc = "multiqc";
#if(! -e "multiqc_fastqc.html"){
#	my $multiqc_cmd = "$multiqc -f ".
#	"-c /export/EC1680U/Pipeline/Trinity/multiqc_fastqc.yaml ".
#	"fastqc_result/ ";
#	system("$multiqc_cmd");
#}

### Trinity de-novo assembly
print CurrTime() . "Start Trinity assembly work...\n";
my $transcripts = "trinity_out_dir/Trinity.fasta";
if(! -e $transcripts){
#	my $Trinity_home = "/usr/local/src/trinityrnaseq-Trinity-v2.4.0";
#	my $docker_basal_cmd = "sudo docker run --rm -v /export/EC1680U:/export/EC1680U -v /export/EC2480U:/export/EC2480U -v $fastq_path:$fastq_path -v `pwd`:`pwd` -w `pwd` kentchendocker/trinityrnaseq:latest";
#	my $trinity_asm_cmd .= "$docker_basal_cmd $Trinity_home/Trinity ".
	my $trinity_asm_cmd .= "Trinity ".
	"--seqType fq --samples_file $input --CPU 32 --max_memory 80G --min_contig_length 300 ".
	"2>&1 | tee $log/Trinity_assembly.log";
	system("$trinity_asm_cmd");
}
print CurrTime() . "Assembly complete!\n\n";
system("sudo chown -R $USER:users trinity_out_dir");

### cd-hit-est for removing redundant assembled transcripts
print CurrTime() . "Start remove redundant transcripts work... (by cd-hit)\n";
if(! -e $transcripts){die "### ERROR: Please make sure assembly is successful.\n";}
my $transcripts_95 = "trinity_out_dir/Trinity.95.fasta";
if(! -e $transcripts_95){
	my $cd_hit_cmd = "cd-hit-est -i $transcripts -o $transcripts_95 ".
#	"-c 0.95 -n 10 -d 0 -M 64000 -T 16 ".
	"-c 0.95 -p 1 -d 0 -b 3 -M 64000 -T 16 ".
	"2>&1 | tee $log/cdhit.log";
	system("$cd_hit_cmd");
}
print CurrTime() . "CD-HIT-EST complete!\n\n";


