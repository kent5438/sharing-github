#! /usr/bin/env perl
use strict;

my $in = (@ARGV > 0) ? shift : 'otu_table.taxID.filter.dense.dedup_L6.txt';
my $out = $in;
$out =~ s/\.txt$//;
$out .= '.charact.txt';

open(OUT,'>'.$out);
open(IN,'<'.$in);
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  my @array = split("\t",$_);

  if ($. == 1) {
    print OUT join("\t",@array)."\n";
    next;
  }

  my $tax = shift(@array);

  $tax =~ s/[^\w-\t\.(?!;$)]+/_/g;

  print OUT join("\t",$tax,@array)."\n";
}
close(IN);
close(OUT);
