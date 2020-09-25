#! /usr/bin/env perl
# Author: kentchen
# Create date: 2017-06-12
# Last update: 2017-10-02

#use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use Time::Seconds;
use lib '/export/EC1680U/kentchen/lib';
use Time qw(CurrTime);

my $pwd = `pwd -L`; chomp($pwd);
my $USER = `id -u -n`;
chomp($USER);

my $tab2xls = "/export/EC1680U/perl/bin/tab2xls.pl"; # need to use for ec2kegg.txt
my $tab2xlsx = "/export/EC1680U/perl/bin/tab2xlsx.pl";

my ($comp, $code);
GetOptions(
	"comp=s"    =>  \$comp,
    "kegg=s"	=> \$code,
);

if (! $comp){
    die "\n### ERROR: The comparison file is necessary\n\n";
} elsif (! $code){
    die "\n### ERROR: KEGG organism code is necessary\n\n";
} else {
    open OUT, ">>CMD_history.txt";
    print OUT CurrTime() . "$0 -c $comp -k $code\n";
    close OUT;
}

### Create folder
my @folders = qw(0_ReadQC 1_AssemblyStats 2_Quantitation 3_DiffExpression 4_Annotation 5_Report);
my ($qc_dir, $stat_dir, $quant_dir, $diff_dir, $anno_dir, $report_dir) = @folders;
my $go_dir = "$anno_dir/GO";
my $kegg_dir = "$anno_dir/KEGG";
my $cog_dir = "$anno_dir/COG";

mkdir("$report_dir");
mkdir("$report_dir/$qc_dir");

### All report variable
my $rsem_report = "$quant_dir/final_RSEM_isoforms_report.txt";
my $de_report = "$diff_dir/final_DE_comparison.report.txt";
my $blastp_report = "$anno_dir/blast_result/uniprot.blastp.report.txt";
my $blastx_report = "$anno_dir/blast_result/uniprot.blastx.report.txt";
my $features_report = "$anno_dir/features/features_report.txt";
my $EC_report = "$kegg_dir/EC_report.txt";
my $COG_report = "$cog_dir/cog_report.txt";
my $GO_report = "$go_dir/Trans_GO.txt";
#my $pathway_report = "$anno_dir/ko2pathway.txt";
my $kegg_report = "$kegg_dir/KEGG_report.txt";


### Merge all reports to one by Rscript
print CurrTime() . "Merged all reports to one xls...\n";

my $report_cmd = 'rsem <- read.table("'.$rsem_report.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'de <- read.table("'.$de_report.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'blastx <- read.table("'.$blastx_report.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="");';
$report_cmd .= 'blastp <- read.table("'.$blastp_report.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="");';
$report_cmd .= 'features <- read.table("'.$features_report.'", header=T, sep="\t", check.names = FALSE, fill = TRUE, quote="\"");';
$report_cmd .= 'EC <- read.table("'.$EC_report.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'EC <- EC[!(EC$`EC number`==""),];';
$report_cmd .= 'COG <- read.table("'.$COG_report.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'GO <- read.table("'.$GO_report.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'KEGG <- read.table("'.$kegg_report.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
#$report_cmd .= 'pathway <- read.table("'.$pathway_report.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'merge1 <- merge(rsem, de, by="transcript_id", all.x=TRUE);';
$report_cmd .= 'merge2 <- merge(merge1, blastx, by="transcript_id", all.x=TRUE);';
$report_cmd .= 'merge3 <- merge(merge2, blastp, by="transcript_id", all.x=TRUE);';
$report_cmd .= 'merge4 <- merge(merge3, features, by="transcript_id", all.x=TRUE);';
$report_cmd .= 'merge5 <- merge(merge4, EC, by="transcript_id", all.x=TRUE);';
$report_cmd .= 'merge6 <- merge(merge5, COG, by="transcript_id", all.x=TRUE);';
$report_cmd .= 'merge7 <- merge(merge6, GO, by="transcript_id", all.x=TRUE);';
#$report_cmd .= 'merge5 <- merge(merge4, pathway, by="KOs", all.x=TRUE);';
$report_cmd .= 'merge8 <- merge(merge7, KEGG, by="transcript_id", all.x=TRUE);';
$report_cmd .= 'write.table(merge8, file="final_report.txt.tmp", sep="\t", quote=F, na="-", row.names = FALSE);';
system("Rscript" . " -e '" . $report_cmd . "'");
#system("awk 'BEGIN{OFS=FS=\"\t\"} {a=\$2; for (i=2; i<NF; i++) \$i=\$(i+1); \$NF=a}1' final_report.txt > final_report.txt.tmp");
system("(head -n 2 final_report.txt.tmp && tail -n +3 final_report.txt.tmp | sort -u -k1,1) > final_report.txt");
system("sed -i 's/\\.x//g' final_report.txt; sed -i 's/\\.y//g' final_report.txt");
system("rm -rf final_report.txt.tmp");

system("$tab2xlsx final_report.txt final_report.xlsx");
system("rm -rf $quant_dir/RSEM_isoforms_length.txt.tmp Trinotate.sqlite");
system("sudo rm -rf tmp.*");

### Transfer specific report files for user
print CurrTime() . "Transfer necessary report files for user...\n";

system("rm -rf 5_Report/4_Annotation/Gene_prediction/work 5_Report/4_Annotation/Gene_prediction/Trinity.95.fasta.transdecoder_dir*");

my $ID = `basename $pwd`; chomp($ID);
my $project_path = "/mnt/NFS/EC2480U-P/scratch/QCResults/fastq/$ID";
if(-e $project_path){
	system("for i in $project_path/reports/multiqc_clean/*; do rsync -av \$i $report_dir/$qc_dir/\$(basename \$i | cut -d'.' -f3,4); done");
} else {
	system("rsync -avq multiqc_data $report_dir/$qc_dir");
	system("rsync -avq multiqc_fastqc.html $report_dir/$qc_dir/multiqc_report.html");
}
# HS19114.clean.multiqc_report_data ---> multiqc_report_data

system("rsync -avq --exclude Trinity.95.fasta.gene_trans_map $stat_dir $report_dir");
system("rsync -avq -m --include='*.xlsx' --include='multiqc*' --include='*mappingStats.txt' --include='bowtie2.sorted_stats/***' --include='*/' --exclude='*' $quant_dir $report_dir");
system("rsync -avq -m --include='*.xlsx' --include='scatter*' --include='*.png' --include='scatter*/***' --include='glimma*/***' --include='*.pdf' --include='*/' --exclude='*' $diff_dir $report_dir");
system("rsync -avq -m --include='*.xlsx' --include='*.pdf' --include='*.png' --include='*.xls' --include='*.txt' --include='*/' --exclude='*' $anno_dir $report_dir");
system("rsync -avq --exclude='Trinity.95.fasta.transdecoder_dir*/' --exclude='work/' $anno_dir/transdecoder/* $report_dir/$anno_dir/Gene_prediction");
#system("rsync -avq $kegg_dir/ec2kegg_README.txt $kegg_dir/Pathway_$code $report_dir/$kegg_dir");
system("rsync -avq final_report.xlsx $report_dir/");
system("rsync -avq /export/EC1680U/kentchen/Template/Help/Transcriptome_Help.pdf $report_dir/Help.pdf");

### Remove long prefix in diff_dir
my @long_prefixs = `find $report_dir/$diff_dir | grep -v "Trans.raw.counts.matrix.xlsx" | grep "Trans.raw.counts.matrix.*"`;
chomp(@long_prefixs);
foreach my $long_prefix (@long_prefixs){
	system("mv $long_prefix \$(echo $long_prefix | sed s/Trans\.raw\.counts\.matrix\.//)");
}


### Generate Report by Rmarkdown

# Copy report_template to current dir)
print CurrTime() . "Copy report_template to current dir...\n";
my $report_template = "/export/EC1680U/Pipeline/Trinity/report_template";
if(! -e "$report_template"){die "\n### ERROR: Please make sure '$report_template' existed!\n\n";}
system("rsync -avq $report_template .");

# Transfer 5_Report to report_template
print CurrTime() . "Transfer 5_Report/* to report_template/files...\n";
system("rsync -avq 5_Report/* report_template/files");

# Build R markdown Script
print CurrTime() . "Build Rmarkdown-Style Report...\n";
chdir("report_template/");
system("perl rmd_forQuant.pl ../comp.txt");
system("perl rmd_forAnno.pl ../comp.txt");
system("Rscript /export/EC1680U/Pipeline/Trinity/build_site.R");
chdir("..");

# Rsync report to root folder & remove report_template
system("rsync -avq report_template/report .");

system("rm -rf report_template report/rmd_forQuant.pl report/rmd_forAnno.pl report/README.txt");

chdir("report/files/4_Annotation/Gene_prediction");
if(-e "Trinity.fasta"){
	system("unlink Trinity.fasta");
	system("ln -s ../../1_AssemblyStats/Trinity.fasta .");
}
if(-e "Trinity.95.fasta"){
	system("unlink Trinity.95.fasta");
	system("ln -s ../../1_AssemblyStats/Trinity.95.fasta .");
}
chdir("$pwd");

print CurrTime() . "Complete!\n\n";

