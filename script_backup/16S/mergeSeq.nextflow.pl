#!/usr/bin/emv perl
use strict;

if ((@ARGV < 1) || (!-e $ARGV[0])) {
  print "Usage: perl $0 mappingFile outputFolder\n";
  exit;
}

my $mapFile = shift;
my $outputFolder = (defined($ARGV[0])?shift(@ARGV):'B.3_demultiplex_quilty_filter');

if (substr($outputFolder,-1) ne '/') {
  $outputFolder .= '/';
}

if (!-e $outputFolder) {
  mkdir($outputFolder);
}

my @sampleOrder;
my %sampleFastqPath;
{
my $sampleNameIdx = 0;
my $filePathIdx = 4;

open(IN,'<'.$mapFile);
chomp($_ = <IN>);
$_ =~ s/^#*//;
my @header = split("\t",$_);
for (my $i=0;$i<@header;$i++) {
  if ($header[$i] eq 'SampleID') {
    $sampleNameIdx = $i;
  }
  elsif ($header[$i] eq 'FilePath') {
    $filePathIdx = $i;
  }
}
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  my @array = split("\t",$_);

  my @filePathArray = split(',',$array[$filePathIdx]);

  push(@sampleOrder,$array[$sampleNameIdx]);
  $sampleFastqPath{$array[$sampleNameIdx]} = \@filePathArray;
}
close(IN);
}

open(OUT,'>'.$outputFolder.'seqs.fna');
foreach my $sampleName (@sampleOrder) {
  my $inputFile = $sampleName.'.extendedFrags.fastq';
  if (!-e $inputFile) {
    #my @fastqPathArray = split(/\//,$sampleFastqPath{$sampleName}[0]);
    $inputFile = $sampleName.'.clean.fastq';
    #$inputFile = $fastqPathArray[$#fastqPathArray];
  }

  if (!-e $inputFile) {
    next;
  }

  my $i = 0;
  open(IN,'<'.$inputFile);
  while (my $name = <IN>) {
    chomp($name);
    chomp(my $seq = <IN>);
    <IN>;
    chomp(my $qual = <IN>);

    $name = substr($name,1);

    print OUT '>'.$sampleName.'_'.$i.' '.$name.' orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0'."\n".$seq."\n";

    $i ++;
  }
  close(IN);
}
close(OUT);
