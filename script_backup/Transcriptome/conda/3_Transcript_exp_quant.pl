#! /usr/bin/env perl
# Author: kentchen
# Create date: 2017-06-12
# Last update: 2017-10-02


use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use Time::Seconds;


my $USER = `id -u -n`;
chomp($USER);
my $rscript = "/usr/bin/Rscript";
my $input;
GetOptions(
    "input=s"  => \$input,
);
if (! $input){
    die "\n### Error: 'input.fofn' is neccessary\n";
} else {
    open OUT, ">>CMD_history.txt";
    print OUT CurrTime() . "$0 -i $input\n";
    close OUT;
}


my $sam_stats = "/export/EC1680U/software/ExpressionAnalysis-ea-utils-27a4809/clipper/sam-stats";
my $tab2xlsx = "/export/EC1680U/perl/bin/tab2xlsx.pl";

### Create folder
my $log = "log";
my @folders = qw(1_AssemblyStats 2_Quantitation 3_DiffExpression);
my ($stat_dir, $quant_dir, $diff_dir) = @folders;
if(! -e $quant_dir){mkdir("$quant_dir");}
if(! -e $log){mkdir("$log");}


### Running Transcript expression quantitation
print CurrTime() . "Start transcript exprssion quantitation work...\n";
my $Trinity_home = "/usr/local/src/trinityrnaseq-Trinity-v2.4.0";
my $docker_basal_cmd = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` trinityrnaseq/trinityrnaseq:latest";
my $trans_95 = "$stat_dir/Trinity.95.fasta";
my $quant_folder_count = `ls $quant_dir | grep -v report | grep -v multiqc | wc -l`; chomp($quant_folder_count);
my $sample_count = `cat $input | wc -l`; chomp($sample_count);

if($quant_folder_count != $sample_count){
	my $trinity_quant_cmd = "$docker_basal_cmd $Trinity_home/util/align_and_estimate_abundance.pl ". 
	"--transcripts $trans_95 --seqType fq --samples_file $input --thread_count 32 ".
	"--est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference ".
	"2>&1 | tee $log/align_and_estimate_abundance.log";
	system("$trinity_quant_cmd");
}

print CurrTime() . "Quantitation Process Complete!\n\n";

### Generate mapping stats & FPKM+1 tmp file & abundance count report
open my $sample_in, "<$input" or die "### ERROR: please check the '$input' existed or not\n";
my @reports;
while(<$sample_in>){
	# mapping stat
	chomp($_);
	my @sp = split/\t/, $_;
	my ($cond, $sample, $r1, $r2) = @sp;
	print CurrTime() . "Generate $sample sam stats report...\n";
	
	my $quant_sample_dir = "$quant_dir/$sample";
	my $sam_stats_cmd = "$sam_stats $quant_sample_dir/bowtie2.bam > $quant_sample_dir/$sample.mappingStats.txt";
	my $tab2xlsx_gene_cmd = "$tab2xlsx $quant_sample_dir/RSEM.genes.results $quant_sample_dir/RSEM.genes.results.xlsx";
	my $tab2xlsx_iso_cmd = "$tab2xlsx $quant_sample_dir/RSEM.isoforms.results $quant_sample_dir/RSEM.isoforms.results.xlsx";
	
	if(! -d "$quant_sample_dir"){
		system("sudo mv -f $sample $quant_dir");
		system("sudo chown -R $USER:$USER $quant_dir");
	}

	system("sudo rm -rf $quant_sample_dir/*.ok");
	if(! -e "$quant_sample_dir/$sample.mappingStats.txt"){system("$sam_stats_cmd");}
	system("$tab2xlsx_gene_cmd");
	system("$tab2xlsx_iso_cmd");

    # Generate sorted.bam & sorted.bam.bai
	my $samtools = "/export/EC1680U/software/samtools-1.3/samtools";
    my $sorted_bam = "$samtools sort -@ 32 $quant_sample_dir/bowtie2.bam -o $quant_sample_dir/bowtie2.sorted.bam";
    system("$sorted_bam");
    my $sorted_bam_index = "$samtools index $quant_sample_dir/bowtie2.sorted.bam";
    system("$sorted_bam_index");


	# Generate FPKM+1 tmp file for logFC calculation
	my $cmd = 'df <- read.table("'.$quant_sample_dir.'/RSEM.isoforms.results", row.names=1, header=T, sep="\t", check.names=F);';
	$cmd .= 'df$FPKM <- df$FPKM + 0.01;';
	$cmd .= 'write.table(df, file="'.$quant_sample_dir.'/RSEM.isoforms.results.plus1.tmp", sep="\t", quote=F, row.names=T, col.names=NA);';
	system($rscript . " -e '" . $cmd . "'");

	# abundance count report
	open my $iso_file, "<$quant_sample_dir/RSEM.isoforms.results";
	open my $iso_report, ">$quant_sample_dir/$sample.RSEM.isoforms.report.txt";
	while(my $iso_txt = <$iso_file>){
		chomp($iso_txt);
		$iso_txt =~ s/expected_count/Count ($sample)/g;
		$iso_txt =~ s/FPKM/FPKM ($sample)/g;
		my @iso_sp = split/\t/, $iso_txt;
		print $iso_report "$iso_sp[0]\t$iso_sp[2]\t$iso_sp[4]\t$iso_sp[6]\n";
	}
	close $iso_file;
	close $iso_report;
	push(@reports, "$quant_sample_dir/$sample.RSEM.isoforms.report.txt");
}
close $sample_in;

### Final abundance count report
my $paste_cmd = "paste @reports > $quant_dir/RSEM.isoforms.report.tmp";
system("$paste_cmd");
my $num_report = scalar @reports;
open my $all_RSEM_report, "<$quant_dir/RSEM.isoforms.report.tmp";
open my $final_RSEM_report, ">$quant_dir/final_RSEM_isoforms_report.txt";
while(<$all_RSEM_report>){
	chomp;
	my @sp = split/\t/, $_;
	my $interval = ($#sp+1)/$num_report;
	print $final_RSEM_report "$sp[0]\t$sp[1]\t";
	for(my $i=2; $i<=$#sp; $i+=$interval){
		print $final_RSEM_report "$sp[$i]\t";
	}
	for(my $j=3; $j<=$#sp; $j+=$interval){
		print $final_RSEM_report "$sp[$j]\t";
	}
	print $final_RSEM_report "\n";
}
close $all_RSEM_report;
close $final_RSEM_report;

system("sudo chown -R $USER:$USER $stat_dir");

system("sed -i 's/\t\$//g' $quant_dir/final_RSEM_isoforms_report.txt");
system("$tab2xlsx $quant_dir/final_RSEM_isoforms_report.txt $quant_dir/final_RSEM_isoforms_report.xlsx");

system("rm -rf $quant_dir/RSEM.isoforms.report.tmp");
system("rm -rf $stat_dir/*bowtie2* $stat_dir/*RSEM*");



sub CurrTime {
    my $t = localtime;
    my $time = "[".$t->hms."]*** ";
    return "$time";
}
