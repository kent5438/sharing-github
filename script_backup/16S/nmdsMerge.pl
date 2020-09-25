#!/usr/bin/env perl
use strict;

if (@ARGV != 2) {
  print 'Usage: perl '.$0.' mappingFile PCoA_weighted_nmds.txt'."\n";
  exit;
}

if ((!-e $ARGV[0]) || (!-e $ARGV[1])) {
  'One of file not exist!'.$/;
  exit;
}

my ($mapping,$nmds) = @ARGV;

my %sampleGroup;
open(IN,'<'.$mapping);
<IN>;
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  my @array = split("\t",$_);
  $sampleGroup{$array[0]} = $array[3];
}
close(IN);

open(OUT,'>PCoA_weighted_nmds.2d.txt');
open(IN,'<'.$nmds);
while (<IN>) {
  chomp;
  if ($_ eq '') {
    last;
  }

  my @array = split("\t",$_);

  if ($array[0] eq "samples") {
    print OUT join("\t",@array,"Description")."\n";
  }
  else {
    print OUT join("\t",@array,$sampleGroup{$array[0]})."\n";
  }
}
close(IN);
close(OUT);
