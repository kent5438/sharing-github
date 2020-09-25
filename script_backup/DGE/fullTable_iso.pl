#!/usr/bin/perl
use strict;

sub uniq {
  my %seen;
  map $seen{$_}=1,@_;
  my @uniq = sort{$a cmp $b} (keys %seen);
  return(@uniq);
}


use Getopt::Long;

my $list = '';
my $outputFolder = '';
my $config = '';

GetOptions(
  'sample:s'  => \$list,
  'config:s' => \$config,
  'output:s' => \$outputFolder
);

if (!$list) {
  $list = 'list.txt';
}

if (!$outputFolder) {
  $outputFolder = 'output';
}

if (!$config) {
  $config = 'dge.cfg';
}

if ((!-e $list) || (!-e $outputFolder) || (!-d $outputFolder) || (!-e $config)) {
  print "Wrong folder!\n";
  exit;
}

my $species = '';
my $mapSpe = 'map';
open(IN,'<'.$config);
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  if ($_ =~ /^genome[=:]/) {
    $species = $';
  }
}
close(IN);

if ($species eq '') {
  print "No species information in config file!\n";
  exit;
}

if (($species eq 'hg19') || ($species eq 'hg38') || ($species eq 'hsa') || ($species eq 'human')) {
  $species = 'human';
  $mapSpe = 'hsa';
}
elsif (($species eq 'mm9') || ($species eq 'mm10') || ($species eq 'mouse')) {
  $species = 'mouse';
  $mapSpe = 'mmu';
}

my %transcriptInfo;
if (($species eq 'human') || ($species eq 'mouse')) {
  #my $file = '/export/EC1680U/daniel/Development/DGE_Version2/mytest/mapTable_'.$species.'.txt';
  my $file = '/export/EC1680U/DataBase/'.(($species eq 'human')?'hg19':'mm10').'/transcriptToGoKegg/mapTable_'.$species.'.txt';

  if (!-e $file) {
    next;
  }

  open(IN,'<'.$file);
  <IN>;
  while (<IN>) {
    chomp;
    if ($_ eq '') {
      next;
    }

    my @array = split("\t",$_);

    $transcriptInfo{$array[0]}{'Description'}{$array[2]} = 1;

    if ($array[3] ne 'NA') {
      $transcriptInfo{$array[0]}{'GO'}{$array[5]}{$array[3]} = 1;
    }

    if ($array[6] ne 'NA') {
      $transcriptInfo{$array[0]}{'Pathway (KEGG)'}{$mapSpe.$array[6]} = 1;
    }

    if ($array[7] ne 'NA') {
      $transcriptInfo{$array[0]}{'EC number'}{$array[7]} = 1;
    }
  }
  close(IN);
}


#my %setInfo;
my @allSample;
open(IN,'<'.$list);
while (<IN>) {
  chomp;
  if (($_ eq '') || ($_ =~ /^#/)) {
    next;
  }

  if ($_ =~ /(SET\d+):(.+)/) {
    my $set = $1;
    my ($test,$ctrl) = split(';',$2);

    my @test_a = split(',',$test);
    my @ctrl_a = split(',',$ctrl);

#    $setInfo{$set}{'test'} = [@test_a];
#    $setInfo{$set}{'ctrl'} = [@ctrl_a];

    push(@allSample,@test_a,@ctrl_a);
  }
}
close(IN);
@allSample = &uniq(@allSample);

my %idToGene;
my %fullTable;
foreach my $sample (@allSample) {
  #open(IN,'<'.$outputFolder.'/1.CalculateExpression/'.$sample.'.genes.results');
  open(IN,'<'.$outputFolder.'/1.CalculateExpression/'.$sample.'.isoforms.results');
  <IN>;
  while (<IN>) {
    chomp;
    if ($_ eq '') {
      next;
    }

    my @array = split("\t",$_);

    my @tids = split(',',$array[0]);

    #$fullTable{$array[0]}{'Effective length'} = $array[3];
    $fullTable{$array[0]}{'Expected count'}{$sample} = $array[4];
    $fullTable{$array[0]}{'TPM'}{$sample} = $array[5];
    $fullTable{$array[0]}{'FPKM'}{$sample} = $array[6];

    foreach my $tid (@tids) {
      $idToGene{'RefSeq'}{$tid} = $array[0];

      foreach my $col ('Description','Pathway (KEGG)','EC number') {
        foreach my $content (keys %{$transcriptInfo{$tid}{$col}}) {
          $fullTable{$array[0]}{$col}{$content} = 1;
        }
      }

      foreach my $ontology ('BP','MF','CC') {
        foreach my $go (keys %{$transcriptInfo{$tid}{'GO'}{$ontology}}) {
          $fullTable{$array[0]}{'GO ('.$ontology.')'}{$go} = 1;
        }
      }
    }
  }
  close(IN);
}

my @allSet;
open(IN,'ls '.$outputFolder.'/3.ebseq/  | ');
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  if ((-e $outputFolder.'/3.ebseq/'.$_) && (-d $outputFolder.'/3.ebseq/'.$_)) {
    push(@allSet,$_);
  }
}
close(IN);

foreach my $set (@allSet) {
  my $file = $outputFolder.'/3.ebseq/'.$set.'/'.$set.'.isoMat.result';
  if (!-e $file) {
    next;
  }

  open(IN,'<'.$file);
  <IN>;
  while (<IN>) {
    chomp;
    if ($_ eq '') {
      next;
    }

    $_ =~ s/\"//g;

    my @array = split("\t",$_);

    #$fullTable{$array[0]}{$set}{'Fold change'} = $array[3];
    #$fullTable{$array[0]}{$set}{'PPEE'} = $array[1];
    $fullTable{$array[0]}{'Fold change'}{$set} = $array[3];
    $fullTable{$array[0]}{'PPEE'}{$set} = $array[1];
  }
  close(IN);
}

open(OUT,'>FullTable.txt');
#print OUT join("\t",'Gene','Effective length');
print OUT join("\t",'Gene');
foreach my $col ('Expected count','TPM','FPKM') {
  foreach my $sample (@allSample) {
    print OUT "\t".$col.' ('.$sample.')';
  }
}
foreach my $col ('Fold change','PPEE') {
  foreach my $set (@allSet) {
    print OUT "\t".$col.' ('.$set.')';
  }
}
print OUT "\t".join("\t",'Description','GO (BP)','GO (MF)','GO (CC)','Pathway (KEGG)','EC number');
print OUT "\n";
foreach my $geneName (sort{$a cmp $b} (keys %fullTable)) {
  #print OUT join("\t",$geneName,$fullTable{$geneName}{'Effective length'});
  print OUT join("\t",$geneName);

  foreach my $col ('Expected count','TPM','FPKM') {
    foreach my $sample (@allSample) {
      my $val = $fullTable{$geneName}{$col}{$sample} + 0;

      print OUT "\t".$val;
    }
  }

  foreach my $col ('Fold change','PPEE') {
    foreach my $set (@allSet) {
      my $val = '--';
      if (defined($fullTable{$geneName}{$col}{$set})) {
        $val = $fullTable{$geneName}{$col}{$set};
      }

      print OUT "\t".$val;
    }
  }

  foreach my $col ('Description','GO (BP)','GO (MF)','GO (CC)','Pathway (KEGG)','EC number') {
    print OUT "\t".join(';',(sort{$a cmp $b} (keys %{$fullTable{$geneName}{$col}})));
  }

  print OUT "\n";
}
close(OUT);

system('perl /export/EC1680U/perl/bin/tab2xlsx.pl FullTable.txt '.$outputFolder.'/5.report/files/FullTable.xlsx');
