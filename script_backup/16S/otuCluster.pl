#!/usr/bin/env perl
use strict;

my $inputFolder = 'B.5_OTU_cluster/97/';
if ((!-e $inputFolder) || (!-d $inputFolder)) {
  print 'Folder '.$inputFolder.' not found!'."\n";
  exit;
}

my %sampleFile;
opendir(my $fh, $inputFolder);
while (readdir($fh)) {
  my $fileName = $_;
  if ($_ =~ /^(.+)\.trimmed\.unique\.opti_mcc\.[\d\.]+\.rep\.fasta/) {
    my $sampleName = $1;
    $sampleFile{$sampleName}{'fasta'} = $fileName;
  }
  elsif ($_ =~ /^(.+)\.trimmed\.unique\.opti_mcc\.[\d\.]+\.cons\.taxonomy/) {
    my $sampleName = $1;
    $sampleFile{$sampleName}{'taxonomy'} = $fileName;
  }
}
closedir($fh);

my $outputFolder = 'report/_reports/data/OTU_cluster/';
if (substr($outputFolder,-1) ne '/') {
  $outputFolder .= '/';
}
if (!-e $outputFolder) {
  system('mkdir -p '.$outputFolder);
}

foreach my $sampleName (keys %sampleFile) {
  my %collect;
  open(IN,'<'.$inputFolder.$sampleFile{$sampleName}{'fasta'});
  while (my $name = <IN>) {
    chomp($name);
    chomp(my $seq = <IN>);

    $name = substr($name,1);

    my @array = split(/[\|\s]+/,$name);

    $collect{$array[1]}{'size'} = $array[2];
    $collect{$array[1]}{'seq'} = $seq;
  }
  close(IN);

  open(IN,'<'.$inputFolder.$sampleFile{$sampleName}{'taxonomy'});
  while (<IN>) {
    chomp;
    if (($_ eq '') || ($. == 1)) {
      next;
    }

    my ($otu,$size,$tax) = split("\t",$_);

    if (defined($collect{$otu}) && defined($collect{$otu}{'size'}) && ($collect{$otu}{'size'} != $size)) {
      print 'Size unequal!'."\n";
      print "\t".$sampleName.' '.$otu."\n";
      exit;
    }

    $collect{$otu}{'tax'} = $tax;
  }
  close(IN);

  open(FA,'>'.$outputFolder.$sampleName.'.otu.fasta');
  open(TXT,'>'.$outputFolder.$sampleName.'.otu.txt');
  print TXT join("\t",'OTU ID','Size','Taxonomy','Sequence')."\n";
  foreach my $otu (sort{$a cmp $b} (keys %collect)) {
    my $size = 0;
    my $seq = '--';
    my $tax = '--';

    if (defined($collect{$otu}{'size'})) {
      $size = $collect{$otu}{'size'};
    }
    if (defined($collect{$otu}{'seq'})) {
      $seq = $collect{$otu}{'seq'};
    }
    if (defined($collect{$otu}{'tax'})) {
      $tax = $collect{$otu}{'tax'};
    }

    print FA '>'.join(' ',$otu,$size,$tax)."\n".$seq."\n";
    print TXT join("\t",$otu,$size,$tax,$seq)."\n";
  }
  close(TXT);
  close(FA);
}
