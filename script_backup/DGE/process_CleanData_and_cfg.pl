#! /usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my ($readdir, $output, $species);
GetOptions(
    "readdir=s" => \$readdir,
    "output=s"  => \$output,
    "species=s"  => \$species,
);

my $pwd = `pwd`; chomp($pwd);
if(!($pwd eq '/export/EC1680U/kentchen/DGE')) {
    die "### ERROR: wrong starting path!\n";
}
if(! -e "$output/reads"){
    system("mkdir -p $output/reads");
}
if(! -e $readdir){
    die "### ERROR: CleanData path is not correct!\n";
}

my @clean_reads = <$readdir/*_001.clean.fastq.gz>;
my $clean_read;
my @first_clean_name=();
my $revised_name;
foreach $clean_read (@clean_reads){
    my ($clean_name, $path, $suffix) = fileparse($clean_read);
#    if($clean_name =~ /(.+)_S[0-9]+_R([0-9])_([0-9]+).clean.fastq.gz/){
    if($clean_name =~ /(.+)_R([0-9])_([0-9]+).clean.fastq.gz/){
        system("rsync -av --ignore-existing $clean_read $output/reads/$1.R$2.fastq.gz");
    }
}
if((! -e "$output/list.txt") or (! -e "$output/dge.cfg")){
	system("rsync -av --ignore-existing list.txt dge.cfg $output");
}
chdir("$output");
my $dge_cfg;
my $dge_line;
open $dge_cfg, "<dge.cfg";
if($species =~ /mouse/){
    while($dge_line = <$dge_cfg>){
        if($dge_line =~ /REFERENCE=\/export\/EC1680U\/DataBase\/(.+)\/RSEM\/reference\/(.+)\.refseq_125polyA/){
		system("sed -i -e 's/$1/mm10/g' dge.cfg");
		system("sed -i -e 's/$2/mm10/g' dge.cfg");
	}
	if($dge_line =~ /NGVEC=\/export\/EC1680U\/DataBase\/(.+)\/RSEM\/reference\/(.+)\.refseq_125polyA\.transcripts\.ngvec/){
		system("sed -i -e 's/$1/mm10/g' dge.cfg");
		system("sed -i -e 's/$2/mm10/g' dge.cfg");
	}
	if($dge_line =~ /genome=(.+)/){
            system("sed -i -e 's/$1/mm10/g' dge.cfg");
        }
    }
} elsif($species =~ /human/){
    while($dge_line = <$dge_cfg>){
	if($dge_line =~ /REFERENCE=\/export\/EC1680U\/DataBase\/(.+)\/RSEM\/reference\/(.+)\.refseq_125polyA/){
                system("sed -i -e 's/$1/hg19/g' dge.cfg");
                system("sed -i -e 's/$2/hg19/g' dge.cfg");
        }
        if($dge_line =~ /NGVEC=\/export\/EC1680U\/DataBase\/(.+)\/RSEM\/reference\/(.+)\.refseq_125polyA\.transcripts\.ngvec/){
                system("sed -i -e 's/$1/hg19/g' dge.cfg");
                system("sed -i -e 's/$2/hg19/g' dge.cfg");
        }
        if($dge_line =~ /genome=(.+)/){
            system("sed -i -e 's/$1/hg19/g' dge.cfg");
        }
    }
}
close $dge_cfg;
system("ln -s $pwd/RNASEQCoreV3_withReport.SR.mod.pl .");

print "\n### Please revise list.txt for the format of 'SET1:<CASE>;<CTRL>'\n
cd to $output, and run  './RNASEQCoreV3_withReport.SR.mod.pl -s list.txt'\n\n";
chdir("reads");
my @revise_reads = <*.fastq.gz>;
my $revise_read;
print "### list.txt\n";
foreach $revise_read (@revise_reads){
	if($revise_read =~ /(.+).R1.fastq.gz/){
		print "$1\n";
	}
}
print "\n\n";
