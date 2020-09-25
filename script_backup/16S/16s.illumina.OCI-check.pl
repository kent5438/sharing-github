#!/usr/bin/env perl
use strict;

chomp(my $pwd = `pwd`);
if (substr($pwd,-1) ne '/') {
  $pwd .= '/';
}

if (!-e $pwd.'data/Mapping.txt') {
  print "Mapping.txt not exist!\n";
  exit;
}

my $cmd = '/usr/local/bin/nextflow run /export/EC1680U/perl/bin/16S/16s.illumina.nf --rm -c /export/EC1680U/perl/bin/16S/16.config';

if ((-e $pwd.'trace.txt') && (-e $pwd.'.nextflow.log')) {
  $cmd .= ' -resume';
}

print $cmd."\n";
system($cmd);

if ($? == 0) {
  print "Exit code: 0\n";
  exit;
}

my @pathFail;
open(IN,'<'.$pwd.'trace.txt');
<IN>;
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  my @array = split("\t",$_);

  my $stat = $array[6];
  my $path = $array[2];

  if ($stat eq 'FAILED') {
    push(@pathFail,$path);
  }
}
close(IN);

my $file = '';
open(IN,'ls work/'.$pathFail[0].'*/slurm.err | ');
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  if (-e $_) {
    $file = $_;
  }
}
close(IN);

my $OCI = 0;
open(IN,'<'.$file);
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  if ($_ =~ /OCI runtime create failed/) {
    $OCI = 1;
  }
}
close(IN);

if ($OCI == 1) {
  $cmd = 'perl '.$0;
  print $cmd."\n";
  exec($cmd);
}
