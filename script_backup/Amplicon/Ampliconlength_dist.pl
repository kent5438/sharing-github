#! /usr/bin/perl

use warnings;
use strict;

my $fastalength = "/export/EC1680U/software/exonerate-2.2.0-x86_64/bin/fastalength";
if(! -e "$ARGV[0].length"){system("$fastalength $ARGV[0] > $ARGV[0].length");}
system('sed -i \'s/\s/\t/g\' '."$ARGV[0].length");

my %hash;
my $over10000count;

open IN, "<$ARGV[0].length";
open OUT, ">$ARGV[0].interval.count.txt";
while(<IN>){
    chomp;
    my @sp = split/\t/, $_;
    my ($length, $contig) = @sp;
    my $index = int($length / 1e3);
    $hash{$index}{count}++;
}
foreach my $index (sort {$a <=> $b} keys %hash) {
    my $count = $hash{$index}{count};
    my $inmil = $index * 1000;
    my $end = $index + 1;
    my $enmil = $end * 1000;
    if($enmil <= 10000){
        print OUT "$inmil-$enmil\t$count\n";
    } else {
        $over10000count += $count;
    }
}
if($over10000count){print OUT "over 10000\t$over10000count\n";}
close IN;
close OUT;

my $pdf = "$ARGV[0].interval.pdf";
my $rscript = "/usr/bin/Rscript";
my $cmd = 'interval <- read.table("'.$ARGV[0].'.interval.count.txt", header=F, sep="\t", check.names=F);';
$cmd .= 'pdf("'.$pdf.'");';
$cmd .= 'barplot(interval$V2, names.arg = interval$V1, las = 2, col = "blue", xlab = "Contig length", ylab = "Count", mgp = c(3, 0.5, 0), cex.lab = 0.7, cex.names = 0.7, cex.axis = 0.7);';
$cmd .= 'dev.off()';
system($rscript . " -e '" . $cmd . "'");
