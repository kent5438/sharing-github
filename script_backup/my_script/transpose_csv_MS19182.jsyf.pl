#!/usr/bin/env perl
use strict;

my @reports;
opendir(my $fh, './');
while (readdir($fh)) {
  if ($_ =~ /^(.+)\.csv/) {
    my $name = $1;

    push(@reports,$name);
  }
}
closedir($fh);

foreach my $name (@reports) {
  my $file = $name.'.csv';
  if (!-e $file) {
    next;
  }

  my %groupSample;
  open(IN,'<'.$file);
  chomp($_ = <IN>);
  my @header = split("\t",$_,-1);
  while (<IN>) {
    chomp;
    if ($_ eq '') {
      next;
    }

    my @array = split("\t",$_,-1);
    for (my $i=0;$i<@header;$i++) {
      if (($array[$i] ne '') && ($_ !~ /^\s+$/)) {
        push(@{$groupSample{$header[$i]}},$array[$i]);
      }
    }
  }
  close(IN);

  open(OUT,'>'.$name.'.txt');
  foreach my $name (sort{$a cmp $b} (keys %groupSample)) {
    print OUT $name.':'.join(',',@{$groupSample{$name}})."\n";
  }
  close(OUT);
}
