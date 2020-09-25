#!/usr/bin/env perl
use strict;

my $pwd = $ENV{'PWD'};
chomp(my $fName = `basename $pwd`);
my $repId = '';
if ($fName =~ /Merge/) {
  $repId = $fName;
  $repId =~ s/^OUTPUT_Merge_//;
}
undef($fName);

my $map = '../data/Mapping.txt';
if ($repId ne '') {
  $map = '../data/Mapping_'.$repId.'.txt';
}
undef($repId);

if (($map eq '') || (!-e $map)) {
  print 'Mapping file '.$map.' not found!'.$/;
  exit;
}

my @samples;
my %sampleName;
open(IN,'<'.$map);
<IN>;
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  my @array = split("\t",$_);

  my $outName = $array[0];
  my ($r1) = split(',',$array[4]);

  chomp(my $fn = `basename $r1`);
  my ($sn) = split(/\.+/,$fn);

  $sampleName{$sn} = $outName;
  push(@samples, $outName);
}
close(IN);

my %stat;
open(IN,'<B.1_QC/multiqc-rawdata/multiqc_data/multiqc_general_stats.txt');
while (<IN>) {
  chomp;
  if (($. == 1) || ($_ eq '')) {
    next;
  }

  my @array = split("\t",$_);

  my $sn = $array[0];
  $sn =~ s/\.R[12]$//;

  my $outName = $sn;
  if (defined($sampleName{$sn})) {
    $outName = $sampleName{$sn};
  }

  $array[4] =~ s/\.0+$//;
  $stat{$outName}{'Clean'} = $array[4];
}
close(IN);

open(IN,'<B.1_QC/multiqc-trimmed/multiqc_data/multiqc_general_stats.txt');
while (<IN>) {
  chomp;
  if (($. == 1) || ($_ eq '')) {
    next;
  }

  my @array = split("\t",$_);

  my $sn = $array[0];
  $sn =~ s/\.R[12]$//;

  $array[4] =~ s/\.0+$//;
  $stat{$sn}{'NoPrimers'} = $array[4];
}
close(IN);

open(IN,'wc -l B.2_join_fastq/*.fastq | ');
while (<IN>) {
  chomp;
  if (($_ eq '') || ($_ !~ /\.fastq$/)) {
    next;
  }

  $_ =~ s/^\s+//;
  $_ =~ s/\s+$//;

  my ($count, $path) = split(/\s+/,$_);

  chomp(my $fn = `basename $path`);
  my ($sn) = split(/\./,$fn);

  my $nSeq = $count/4;

  $stat{$sn}{'Joined'} = $nSeq;
}
close(IN);

opendir(my $fh, 'B.4_removing_chimeras/');
while (readdir($fh)) {
  if ($_ !~ /\.trimmed\.fa/) {
    next;
  }

  my $fn = $_;
  my ($sn) = split(/\./,$fn);

  chomp(my $nl = `grep -c "^>" B.4_removing_chimeras/${fn}`);
  $nl =~ /(\d+)/;
  $nl = $1;

  $stat{$sn}{'NoChimeras'} = $nl;
}
closedir($fh);

=h
open(IN,'grep -c "^>" B.4_removing_chimeras/*.trimmed.fa | ');
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  $_ =~ s/^\s+//;
  $_ =~ s/\s+$//;

  my ($path, $count) = split(':',$_);

  chomp(my $fn = `basename $path`);
  my ($sn) = split(/\./,$fn);

  $stat{$sn}{'NoChimeras'} = $count;
}
close(IN);
=cut

foreach my $sn (keys %stat) {
  my $file = '../OUTPUT/B_statistic_analysis/'.$sn.'.statistic_sample.txt';
  if (!-e $file) {
    next;
  }

  my @header;
  open(IN,'<'.$file);
  while (<IN>) {
    chomp;
    if ($_ eq '') {
      next;
    }

    if ($. == 1) {
      @header = split("\t",$_);
      next;
    }

    my @array = split("\t",$_);

    for (my $i=1;$i<@array;$i++) {
      $stat{$sn}{$header[$i]} = $array[$i];
    }
  }
  close(IN);
}

my @cols = qw / Clean   NoPrimers       Joined  NoChimeras      Total_Tags      Unique_Tags     Taxon_Tags      OTUs / ;

open(OUT,'>B_statistic_analysis/reads-tags_summary.txt');
print OUT join("\t",'SampleID',@cols)."\n";
foreach my $sn (@samples) {
  print OUT $sn;

  foreach my $col (@cols) {
    my $val = '--';
    if (defined($stat{$sn}{$col})) {
      $val = $stat{$sn}{$col};
    }

    print OUT "\t".$val;
  }

  print OUT "\n";
}
close(OUT);
