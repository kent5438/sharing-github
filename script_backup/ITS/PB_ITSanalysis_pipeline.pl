#! /usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use Time::Seconds;

my $t = localtime;
my $time = "[".$t->hms."] ";

### DB
my $ref2ChimeraUchime="/export/EC1680U/metagenomics/ITSdb/chimera/uchime_reference_dataset_28.06.2017/uchime_reference_dataset_28.06.2017.fasta";
if (! -e "$ref2ChimeraUchime"){die "\n### ERROR: $ref2ChimeraUchime not found!!\n\n"}

my $refAlignDB_97="/export/EC1680U/metagenomics/ITSdb/align/UNITEv6_sh_97.aln.fasta";
my $refAlignDB_99="/export/EC1680U/metagenomics/ITSdb/align/UNITEv6_sh_99.aln.fasta";
#if (! -e "$refAlignDB_99"){die "\n### ERROR: $refAlignDB_99 not found!!\n\n"}

my $refFastaDB_97="/export/EC1680U/metagenomics/ITSdb/mothur/UNITEv6_sh_97.fasta";
my $refFastaDB_99="/export/EC1680U/metagenomics/ITSdb/mothur/UNITEv6_sh_99.fasta";
if (!-e "$refFastaDB_97"){die "### ERROR: $refFastaDB_97 not found!!\n\n"}

my $refTaxonomy_97="/export/EC1680U/metagenomics/ITSdb/mothur/UNITEv6_sh_97.tax"; 
my $refTaxonomy_99="/export/EC1680U/metagenomics/ITSdb/mothur/UNITEv6_sh_99.tax"; 
if (!-e "$refTaxonomy_97"){die "### ERROR: $refTaxonomy_97 not found!!\n\n"}


### Create folders
my @folders = qw(0_cleanFastq 1_joinPairs 2_removePrimers 3_screenQuality 4_checkChimera 5_assignTaxa 6_biom metaITS_OUTPUT);
my ($clean_dir, $join_dir, $rmPrimer_dir, $quality_dir, $chimera_dir, $taxa_dir, $biom_dir, $meta_dir) = @folders;
my $eachOTU_dir = "$meta_dir/OTUs_eachSample";
foreach my $folder (@folders){
	if(! -e $folder){mkdir "$folder";}
}


### Setting options
my $USER = `id -u -n`;
chomp($USER);
my ($help, $input, $primer, $sequencer);
GetOptions(
	"help!"	=> \$help,
	"input=s"	=> \$input,
	"primer=s"	=> \$primer,
	"sequencer=s"	=> \$sequencer,
);

if($help){
	die "\n### USAGE ###
* MiSeq / Pacbio ITS pipeline
[necessary files in project root]:
- input.fofn				(read file)
- ITS1-2primers.alienSeq	(primer)
- PB17028_barcode.txt		(barcode)
- report-ITS.cfg			(report config)

e.g.)
% PB_ITSanalysis_pipeline.pl -i input.fofn -p primer.fa -s PB

-i|--input:
Sample information & Clean read repository path
e.g (need to be fastq format, instead of fastq.gz)
condA	sample_1	0_cleanFastq/sample_1.R1.fastq	0_cleanFastq/sample_1.R2.fastq
condB	sample_2	0_cleanFastq/sample_2.R1.fastq	0_cleanFastq/sample_2.R2.fastq
condC	sample_3	0_cleanFastq/sample_3.R1.fastq	0_cleanFastq/sample_3.R2.fastq

-p|--primer: 'ITS1-2primers.alienSeq' is necessary which need to contain ITS primer
-s|--sequencer: PB / MS\n
";
} elsif (! $input){
    die "\n### Error: 'input.fofn' is neccessary\n";
} elsif (! $primer){
	die "\n### Error: 'ITS1-2primers.alienSeq' is neccessary\n";
} elsif (! (($sequencer =~ /PB/) or ($sequencer =~ /MS/))){
	die "\n### Error: -s|--sequencer option needs to be 'MS|PB'\n";
}

print "\n[Reminder]: please set your 'report-ITS.cfg' correctly!\n\n";
sleep(5);

### Try to detect blank line, or the count of samples will be false
open IN, "<$input" or die "### ERROR: please check the '$input' existed or not\n";
while(<IN>){
	chomp;
	if(/^\s*$/){die "\n### ERROR: $input detected blank line!\n\n";}
}
close IN;

my $sample_count = `cat $input | wc -l`;
chomp($sample_count);

### Check reads are matched with input.fofn & Generate Sample_list.txt
open my $sample_in, "<$input" or die "### ERROR: please check the '$input' existed or not\n";
open my $sample_list, ">Sample_list.txt";

my $log = "log";
if(! -e $log){mkdir("$log");}

### Main pipeline
my $pwd = `pwd -L`;
chomp($pwd);

while(<$sample_in>){
	chomp($_);
	my @sp = split/\t/, $_;
	my ($cond, $sample, $r1, $r2) = @sp;
	if(! -f "$r1"){die "\n### ERROR: Cannot find $r1!\nYour read files need to be concordant with $input\n";}
	print $sample_list "$sample\n";
	
	my $dash_sample = $sample;
	if($dash_sample =~ /-/){$dash_sample =~ s/-/\\-/g;} # prevent dash bug by mothur

	if($sequencer =~ /MS/){
		print "$time* Detect this case is MiSeq paired-end reads, starting $sample PEAR merged work...\n";
		if (! $r2) {die "\n### ERROR: Please make sure Read_2 are existed!\n";}
## Merge paired-end reads by PEAR if paired-end data
		my $pear_cmd = "/export/EC1680U/software/anaconda2/bin/pear -m 600 -n 50 -s 2 -u 1 -p 0.05 -v 15 -j 24 ".
		"-f $r1 -r $r2 -o $join_dir/$sample 2>&1 | tee $join_dir/$sample.pear.log";
		if(! -e "$join_dir/$sample\.assembled.fastq"){
			system("$pear_cmd");
		}
	} elsif($sequencer =~ /PB/){
		print "$time* Detect this case is Pacbio Single-end reads, skip $sample PEAR merged work...\n";
		if ($2) {die "\n### ERROR: Pacbio must be single end data!\n or is there containing any space?\n";}
		if(! -e "$join_dir/$sample\.assembled.fastq"){
			system("ln -s $pwd/$r1 $pwd/$join_dir/$sample\.assembled.fastq");
		}
	} else {
		die "### ERROR: please check $r1 & $r2 are located in correct directory!\n";
	}
	print $sample_list "$sample\n";

## Primer trimming by AlienTrimmer
	print "$time* $sample\'s Primer trimming work by AlienTrimmer...\n";
	my $trim_cmd = "/export/EC1680U/software/AlienTrimmer_0.4.0/src/AlienTrimmer -q 0 -l 300 -p 50 ".
	"-i $join_dir/$sample\.assembled.fastq -o $rmPrimer_dir/$sample\.assembledTrimPrimers.fastq ".
	"-c $primer 2>&1 | tee $rmPrimer_dir/$sample\.assembledTrimPrimers.log";
	if(! -e "$rmPrimer_dir/$sample\.assembledTrimPrimers.fastq"){
		system("$trim_cmd");
	}

## Quality improvement by Trimmomatic
	print "$time* $sample\'s Quality improvement work by Trimmomatic...\n";
	my $trimmomatic_cmd = "/export/EC1680U/software/anaconda2/bin/trimmomatic SE -threads 24 -phred33 ".
	"$rmPrimer_dir/$sample\.assembledTrimPrimers.fastq $quality_dir/$sample\.qualified.fastq ".
	"SLIDINGWINDOW:10:10 MINLEN:200 2>&1 | tee $quality_dir/$sample\.trimmomatic.log";
	if(! -e "$quality_dir/$sample\.qualified.fastq"){
		system("$trimmomatic_cmd");
	}
## Quality info by Mothur
	print "$time* $sample\'s Quality info work by Mothur...\n";
	if(! -e "$quality_dir/$sample.qualified.fasta"){
		if($sequencer =~ /PB/){
			my $mothur_info_PB = "/export/EC1680U/software/anaconda2/bin/mothur \"#fastq.info(pacbio,fastq=$quality_dir/$sample.qualified.fastq)\"";
			system("$mothur_info_PB");
		}
		if($sequencer =~ /MS/){
			my $mothur_info_MS = "/export/EC1680U/software/anaconda2/bin/mothur \"#fastq.info(format=illumina1.8+,qfile=T,fastq=$quality_dir/$sample.qualified.fastq)\"";
			system("$mothur_info_MS");
		}
	}
## Remove Chimera by Mothur (beware mothur chimera.uchime will be crashed if sample name contains dash)
	print "$time* $sample\'s Remove Chimera work by Mothur...\n";
	my $mothur_uchime = "/export/EC1680U/software/anaconda2/bin/mothur ".
	"\"#chimera.uchime(processors=32, fasta=$quality_dir/$dash_sample.qualified.fasta, reference=$ref2ChimeraUchime)\" 2>&1 | tee $chimera_dir/$sample.chimera.log";
	my $mothur_rmSeq = "/export/EC1680U/software/anaconda2/bin/mothur ".
	"\"#remove.seqs(accnos=$quality_dir/$sample.qualified.ref.uchime.accnos, fasta=$quality_dir/$sample.qualified.fasta, qfile=$quality_dir/$sample.qualified.qual)\"";
	if(! -e "$quality_dir/$sample.qualified.ref.uchime.chimeras"){
		system("$mothur_uchime");
	}
	if(! -e "$chimera_dir/$sample.effective.fasta"){
		system("$mothur_rmSeq");
		system("mv $quality_dir/$sample.qualified.pick.fasta $chimera_dir/$sample.effective.fasta");
		system("mv $quality_dir/$sample.qualified.pick.qual $chimera_dir/$sample.effective.qual");
	}

# fasta & qual -> fastq to calculate stats easily, discard fastq lately?
	my $fastaQual2fastq = "/export/EC1680U/software/anaconda2/bin/mothur \"#make.fastq(fasta=$chimera_dir/$sample.effective.fasta, qfile=$chimera_dir/$sample.effective.qual)\"";
	if(! -e "$chimera_dir/$sample.effective.fastq"){
		system("$fastaQual2fastq");
	}
	system("rm -rf $quality_dir/$sample.*.tmp.uchime_formatted $quality_dir/$sample.*.uchime*.temp $quality_dir/$sample*.tmp");
## Statistics info from above processing data
	print "$time* $sample\'s Statistics info from above processing data...\n\n";


}
close $sample_in;
close $sample_list;




























