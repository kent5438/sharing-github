#!/usr/bin/perl
#author: JSYF
#modified date: 2018/07/03
use strict;

sub generateFullTable {
  my @allSample = @{(shift)};
  my %species;
  $species{'Species'} = shift;
  my %folders;
  #$folders{'Output'} = shift;
  #$folders{'CalculateExpression'} = $folders{'Output'}.'/1.CalculateExpression/';
  #$folders{'Ebseq'} = $folders{'Output'}.'/3.ebseq/';
  #$folders{'StatsOut'} = $folders{'Output'}.'/4.stats_out/';
  #$folders{'Report'} = $folders{'Output'}.'/5.report/';
  $folders{'CalculateExpression'} = shift;
  $folders{'Ebseq'} = shift;
  $folders{'StatsOut'} = shift;

  if (($species{'Species'} eq 'hg19') || ($species{'Species'} eq 'hg38') || ($species{'Species'} eq 'hsa') || ($species{'Species'} eq 'human')) {
    $species{'Species'} = 'human';
    $species{'Spe'} = 'hsa';
    $species{'Genome'} = 'hg19';
  }
  elsif (($species{'Species'} eq 'mm9') || ($species{'Species'} eq 'mm10') || ($species{'Species'} eq 'mmu') || ($species{'Species'} eq 'mouse')) {
    $species{'Species'} = 'mouse';
    $species{'Spe'} = 'mmu';
    $species{'Genome'} = 'mm10';
  }
  elsif (($species{'Species'} eq 'fruitfly') || ($species{'Species'} eq 'dm6') || ($species{'Species'} eq 'dme')) {
  	$species{'Species'} = 'fruitfly';
	$species{'Spe'} = 'dme';
    $species{'Genome'} = 'dm6';
  }
  else {
    $species{'Spe'} = 'map';
    $species{'Genome'} = $species{'Species'};
  }

  my %files;
#  $files{'TranscriptInfo'} = '/export/EC1680U/DataBase/'.$species{'Genome'}.'/transcriptToGoKegg/mapTable_'.$species{'Species'}.'.txt';
  $files{'TranscriptInfo'} = '/export/EC1680U/DataBase/D_melanogaster/mapTable_Dm.txt';

  foreach my $folder (keys %folders) {
    $folders{$folder} .= '/';
    $folders{$folder} =~ s/\/+/\//g;

    if (!-e $folders{$folder}) {
      system('mkdir -p '.$folders{$folder});
    }
  }

  my %transcriptInfo;
  if (-e $files{'TranscriptInfo'}) {
    open(IN,'<'.$files{'TranscriptInfo'});
    while (<IN>) {
      chomp;
      if ($_ eq '') { next; }

      my @array = split("\t",$_);

      if ($array[2] ne 'NA') {
        $transcriptInfo{$array[0]}{'Description'}{$array[2]} = 1;
      }
      if ($array[3] ne 'NA') {
        $transcriptInfo{$array[0]}{'GO'}{$array[5]}{$array[3]} = 1;
      }

      if ($array[6] ne 'NA') {
        $transcriptInfo{$array[0]}{'Pathway (KEGG)'}{$species{'Spe'}.$array[6]} = 1;
      }

      if ($array[7] ne 'NA') {
        $transcriptInfo{$array[0]}{'EC number'}{$array[7]} = 1;
      }
    }
    close(IN);
  }

  my %fullTable;
  foreach my $sample (@allSample) {
    open(IN,'<'.$folders{'CalculateExpression'}.$sample.'.genes.results') or (die($folders{'CalculateExpression'}.$sample.'.genes.results'." open failed\n"));
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
        #$idToGene{'RefSeq'}{$tid} = $array[0];

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
  open(IN,' ls '.$folders{'Ebseq'}.' |');
  while (<IN>) {
    chomp;
    if ($_ eq '') {
      next;
    }

    if ((-e $folders{'Ebseq'}.$_) && (-d $folders{'Ebseq'}.$_)) {
      push(@allSet,$_);
    }
  }
  close(IN);

  foreach my $set (@allSet) {
    my $file = $folders{'Ebseq'}.$set.'/'.$set.'.GeneMat.results';
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

      $fullTable{$array[0]}{'Fold change'}{$set} = $array[3];
      $fullTable{$array[0]}{'PPEE'}{$set} = $array[1];
    }
    close(IN);
  }

  open(OUT,'>'.$folders{'StatsOut'}.'FullTable.txt');
  print OUT 'Gene';
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
}

if (@ARGV < 2) {
  print "Usage: perl $0 outputFolder DGE.json\n";
  exit;
}

my ($folder,$json) = @ARGV;
if (substr($folder,-1) ne '/') {
  $folder .= '/';
}
$folder .= 'results/';

if ((!-e $folder.'1.CalculateExpression/') || (!-e $folder.'3.ebseq/') || (!-e $folder.'4.stats_out/')) {
  print "Wrong folder!\n";
  exit;
}
elsif (!-e $json) {
  print "JSON file $json not found!\n";
  exit;
}

my $genome;
my @samples;
open(IN,'<'.$json);
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  if ($_ =~ /^\s*\"genome\":\s*\"([^\"]+)\"/) {
    $genome = $1;
  }
  elsif ($_ =~ /^\s*\"Sample_Name\":\s*\"([^\"]+)\"/) {
    my $sampleName = $1;
    push(@samples,$sampleName);
  }
}
close(IN);

if (($genome eq '') || (@samples == 0)) {
  print "No genome or sample information in json file!\n";
  exit;
}

&generateFullTable(\@samples,$genome,$folder.'1.CalculateExpression/',$folder.'3.ebseq/',$folder.'4.stats_out/');

