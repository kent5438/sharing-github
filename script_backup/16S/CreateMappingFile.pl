#! /usr/bin/perl

use strict;
use warnings;

my ($readType, $barcode) = @ARGV;
if (@ARGV == 0){die "Usage: $0 <PE|SE> <PS17015-16S-barcode_1.txt>\n";}
if(!($readType =~ /PE|SE/)){die "### ERROR: 
Usage: $0 <PE|SE> <PS17015-16S-barcode_1.txt>\n";}
elsif(! $barcode){print "### [Warnings]: barcode list is used for rename Pacbio CCS file name\n";}

my $pwd = `pwd -L`; chomp($pwd);
if(! ($pwd =~ /\/data/)){die "### ERROR: Please make sure you are located at /path/of/project/data\n";}

my $checkFastq = `ls *.fastq | wc -l`; chomp($checkFastq);
if($checkFastq == 0){die "### ERROR: No Fastq file in data/ !\n";}


open OUT, ">Mapping.txt";

print OUT "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\tFilePath\n";

if($readType eq "PE"){
	my @fastqs = `ls *.R1.clean.fastq | grep -v Undetermined | sort -V -u`; chomp(@fastqs);
	foreach my $fastq (@fastqs){
		my @sp = split/\.R1/, $fastq;
		my $prefix = $sp[0];
		(my $noDash_prefix = $prefix) =~ s/-/_/g;
		my $R1_fastq = "$prefix\.R1.clean.fastq";
		my $R2_fastq = "$prefix\.R2.clean.fastq";
		print OUT "$noDash_prefix\t\t\t$prefix\t$pwd/$R1_fastq,$pwd/$R2_fastq\n";
	}
	system("cp /export/EC1680U/Pipeline/16S/MiSeq-16S_multiqc_config.yaml ../");
#	system("cp /export/EC1680U/Pipeline/16S/V3V4_primer.fa primer.fa");
	system("cp /export/EC1680U/perl/bin/16S/16s.illumina.v4.nf ../");
	system("cp /export/EC1680U/perl/bin/16S/16.v4.config ../");
	print "\n### [Notice]: Beware of the 4th Group column!
You could modify them if necessary...\n\n";
	
	system("perl /export/EC1680U/perl/bin/16S/subsample_fastq.pl");
}

if($readType eq "SE"){
#	system("perl /export/EC1680U/perl/bin/16S/PBccs_renameForNew16S.pl $barcode");
	
	# need to reformat phred quality to avoid the bug of qiime split_library
	mkdir "ccs.raw";
	system("for i in `ls *.fastq | awk -F'.fastq' '{print \$1}'`; do \
		/export/EC1680U/software/bbmap/reformat.sh in=\${i}.fastq out=\${i}.fixed.fastq ; \
		mv \${i}.fastq ccs.raw; \
		mv \${i}.fixed.fastq \${i}.fastq; \
		done");

	my @fastqs = `ls *.fastq | sort -V -u`; chomp(@fastqs);
	foreach my $fastq (@fastqs){
		my @sp = split/\./, $fastq;
		my $prefix = $sp[0];
		(my $noDash_prefix = $prefix) =~ s/-/_/g;
		print OUT "$noDash_prefix\t\t\t$prefix\t$pwd/$fastq\n";
	}
	system("cp /export/EC1680U/Pipeline/16S/PB-16S_multiqc_config.yaml ../");

	system("cp /export/EC1680U/perl/bin/16S/16s.PacBio.v2.nf ../");
#	system("cp /export/EC1680U/Pipeline/16S/PB_primer.fa primer.fa");
	system("cp /export/EC1680U/perl/bin/16S/16.v4.config ../");
	print "\n### [Notice]: Beware of the 4th Group column!
You could modify them if necessary...\n\n";
}
close OUT;

#system("cp /export/EC1680U/perl/bin/16S/data_preprocess.pl ../");
