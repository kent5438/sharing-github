#!/usr/bin/env perl
use strict;

my $table = shift;

my %rawName;
open(IN,'</references/SILVA_132_release/taxonomy/16S_only/99/raw_taxonomy.txt.bak');
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  my ($taxId, $tax) = split("\t",$_);

  my @taxonomy = split(';',$tax);

  $rawName{$taxId} = \@taxonomy;
}
close(IN);

open(OUT,'>motu_table.taxID.filter.dense.dedup.txt');
open(IN,'<'.$table);
while (<IN>) {
  chomp;
  if ($_ =~ /^#/) {
    print OUT $_."\n";
    next;
  }

  my @array = split("\t",$_);
  my $taxId = shift(@array);
  my $tax = pop(@array);

  my @taxonomy = split('; ',$tax);

  my @ref = @taxonomy;
  if (defined($rawName{$taxId})) {
    @ref = @{$rawName{$taxId}};
  }

  for (my $i=0;$i<@taxonomy;$i++) {
    if ($taxonomy[$i] !~ /_unclassified$/) {
      # replace by rawName
      $taxonomy[$i] = $ref[$i];
    }

    $taxonomy[$i] =~ s/^[kpcofgs]__//;
  }

  print OUT join("\t", $taxId, @array, join('; ',@taxonomy))."\n";
}
close(IN);
close(OUT);
