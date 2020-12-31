#! /usr/bin/perl

use strict;

open IN, "<Mapping.txt";
while(<IN>){
    chomp($_);
    my @sp = split/\t/, $_;
    my ($cond, $sample, $r1, $r2) = @sp;
    if(! -f "$r1"){die "\n### ERROR: Cannot find $r1!\n";}

    my $dash_sample = $sample;
        if (! $r2) {die "\n### ERROR: Please make sure Read_2 are existed!\n";}

## Merge paired-end reads by PEAR if paired-end data
        my $pear_cmd = "/export/EC1680U/software/anaconda2/bin/pear -m 400 -n 50 -s 2 -u 1 -p 0.05 -v 15 -j 24 ".
        "-f $r1 -r $r2 -o $sample 2>&1 | tee $sample.pear.log";
        if(! -e "$sample\.assembled.fastq"){
            system("$pear_cmd");
        }
}
close IN;

system("cat *.pear.log | grep 'Assembled reads \\.' | awk -F': ' '{print \$2}' | awk -F' / ' '{print \$1}' > count.txt");
system("cat *.pear.log | grep 'Assembled reads file\\.' | awk -F': ' '{print \$2}' | awk -F' / ' '{print \$1}' | cut -d'.' -f1 | sed 's/-//g' > sample.txt");
system("paste sample.txt count.txt | sort -n -k1,1 > all.count.txt");

