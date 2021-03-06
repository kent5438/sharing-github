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

### Create folder
my $log = "log";
my @folders = qw(1_AssemblyStats 2_Quantitation 3_DiffExpression);
my ($stat_dir, $quant_dir, $diff_dir) = @folders;
if(! -e $stat_dir){mkdir("$stat_dir");}
if(! -e $log){mkdir("$log");}
### Transfer completed assembled transcripts file
my $USER = `id -u -n`;
chomp($USER);
system("cp trinity_out_dir/Trinity.fasta trinity_out_dir/Trinity.95.fasta $stat_dir");

### Transcript length generate
chdir("$stat_dir");
system("fastalength Trinity.95.fasta > Trinity.95.length.txt");
system('sed -i \'s/\s/\t/g\' Trinity.95.length.txt');

### Generate transcripts interval distribution
print CurrTime() . "Generate transcripts interval distribution...\n";
my %hash;
my $over3000count;

open IN, "<Trinity.95.length.txt";
open OUT, ">Trinity.95.interval.count.txt";
while(<IN>){
	chomp;
	my @sp = split/\t/, $_;
	my ($length, $contig) = @sp;
	my $index = int($length / 1e2);
	$hash{$index}{count}++;
}
foreach my $index (sort {$a <=> $b} keys %hash) {
	my $count = $hash{$index}{count};
	my $inmil = $index * 100;
	my $end = $index + 1;
	my $enmil = $end * 100;
	if($inmil < 3000){
		print OUT "$inmil-$enmil\t$count\n";
	} else {
		$over3000count += $count;
	}
}
if($over3000count){print OUT "over 3000\t$over3000count\n";}
close IN;
close OUT;

### Plot by R
print CurrTime() . "Plot transcripts interval length distribution...\n";
my $cmd = 'interval <- read.table("Trinity.95.interval.count.txt", header=F, sep="\t", check.names=F);';
$cmd .= 'png("Trinity.95.dist.png", width=10, height=10, units="in", res=300);';
$cmd .= 'barplot(interval$V2, names.arg = interval$V1, las = 2, col = "blue", xlab = "Transcript length", ylab = "Count", mgp = c(3, 0.5, 0), cex.lab = 0.7, cex.names = 0.7, cex.axis = 0.7);';
$cmd .= 'dev.off()';
system("Rscript" . " -e '" . $cmd . "'");


### Trinity assembly evaluation
print CurrTime() . "Trinity assembly evaluation...\n";
my $transcript_stats_cmd = "AssemblyStatistics contigs Trinity.95.fasta";
system("$transcript_stats_cmd");
system("grep ': ' Trinity.95.fastaStats.txt | sed 's/: /\t/g' > Trinity_assembled.final.stats.txt");
system("grep -v ': ' Trinity.95.fastaStats.txt | sed 's/ /\t/g' >> Trinity_assembled.final.stats.txt");
system("sed -i '1i ### skip this line' Trinity_assembled.final.stats.txt");

chdir("..");


