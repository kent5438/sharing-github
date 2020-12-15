#!/usr/bin/env perl
use strict;

sub uniq {
  my %seen;
  map $seen{$_}=1,@_;
  my @uniq = sort{$a cmp $b} (keys %seen);
  return(@uniq);
}

sub readBiom {
  use JSON;

  my $input = shift;

  if (!-e $input) {
    &usage();
    exit;
  }

  my $strJson = '';
  open(IN,'<'.$input);
  while (<IN>) {
    chomp;
    if ($_ eq '') {
      next;
    }

    $strJson .= $_;
  }
  close(IN);

  my $refJson = decode_json($strJson);
  return($refJson);
}

sub writeBiom {
  use JSON;

  my ($ref,$output) = @_;

  my $strJson = to_json($ref,{pretty => 1});

  open(OUT,'>'.$output);
  print OUT $strJson."\n";
  close(OUT);

  return(0);
}

## check biom
sub checkBiom { # arg: inputFileName
  my $input = shift;

  my $jsonHashRef = &readBiom($input);

  my $nRow = scalar(@{$jsonHashRef->{'rows'}});
  my $nCol = scalar(@{$jsonHashRef->{'columns'}});

  my $nTableRow = scalar(@{$jsonHashRef->{'data'}});
  my $nTableCol0 = scalar(@{$jsonHashRef->{'data'}[0]});
  my $nTableCol1 = scalar(@{$jsonHashRef->{'data'}[1]});
  my $nTableCol2 = scalar(@{$jsonHashRef->{'data'}[2]});

  print join("\n",$nRow,$nCol,$nTableRow,$nTableCol0,$nTableCol1,$nTableCol2)."\n";
}

## task: uniqByTax
sub uniqByTax { # arg: inputFileName, outputFileName
  my ($input,$output) = @_;

  if (scalar (@{[split(',',$input)]}) > 1) {
    print "Invalid number of input file\n";
    &usage();
    exit;
  }
  if (scalar (@{[split(',',$output)]}) > 1) {
    print "Invalid number of input file\n";
    &usage();
    exit;
  }

  my $jsonHashRef = &readBiom($input);
  my %taxIdIdx;
  for (my $i=0;$i<$jsonHashRef->{'shape'}[0];$i++) {
    if (!defined($jsonHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'})) {
      next;
    }
    my $tax = '';
    if (ref($jsonHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'}) eq 'ARRAY') {
      $tax = join(';;',@{$jsonHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'}});
    }
    else {
      $tax = $jsonHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'};
    }

    push(@{$taxIdIdx{$tax}{'idx'}},$i);

    if (!defined($taxIdIdx{$tax}{'id'}) || ($taxIdIdx{$tax}{'id'} eq '')) {
      $taxIdIdx{$tax}{'id'} = $jsonHashRef->{'rows'}[$i]{'id'};
    }
  }

  my @taxOrder = sort{$taxIdIdx{$a}{'id'} cmp $taxIdIdx{$b}{'id'}} (keys %taxIdIdx);

  my @newRows;
  foreach my $tax (@taxOrder) {
    my %content;
    $content{'id'} = $taxIdIdx{$tax}{'id'};
    @{$content{'metadata'}{'taxonomy'}} = split(';;',$tax);

    push(@newRows,\%content);
  }

  my @newData;
  for (my $iNew=0;$iNew<@taxOrder;$iNew++) {
    my $tax = $taxOrder[$iNew];
    foreach my $i (@{$taxIdIdx{$tax}{'idx'}}) {
      for (my $j=0;$j<@{$jsonHashRef->{'data'}[$i]};$j++) {
        $newData[$iNew][$j] += $jsonHashRef->{'data'}[$i][$j];
      }
    }
  }

  my %outputHash = %{$jsonHashRef};
  $outputHash{'rows'} = \@newRows;
  $outputHash{'data'} = \@newData;
  $outputHash{'shape'}[0] = scalar @newRows;

  &writeBiom(\%outputHash,$output);
}

## task: uniqById
sub uniqById { # arg: inputFileName, outputFileName
  my ($input,$output) = @_;

  if (scalar (@{[split(',',$input)]}) > 1) {
    print "Invalid number of input file\n";
    &usage();
    exit;
  }
  if (scalar (@{[split(',',$output)]}) > 1) {
    print "Invalid number of input file\n";
    &usage();
    exit;
  }

  my $jsonHashRef = &readBiom($input);

  my %taxidTaxIdx;
  for (my $i=0;$i<$jsonHashRef->{'shape'}[0];$i++) {
    if (!defined($jsonHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'})) {
      next;
    }
    my $tax = '';
    if (ref($jsonHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'}) eq 'ARRAY') {
      $tax = join(';;',@{$jsonHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'}});
    }
    else {
      $tax = $jsonHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'};
    }

    push(@{$taxidTaxIdx{$jsonHashRef->{'rows'}[$i]{'id'}}{'idx'}},$i);

    if (!defined($taxidTaxIdx{$jsonHashRef->{'rows'}[$i]{'id'}}{'tax'}) || ($taxidTaxIdx{$jsonHashRef->{'rows'}[$i]{'id'}}{'tax'} eq '')) {
      $taxidTaxIdx{$jsonHashRef->{'rows'}[$i]{'id'}}{'tax'} = $jsonHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'};
    }
  }

  my @taxidOrder = sort{$a cmp $b} (keys %taxidTaxIdx);

  my @newRows;
  foreach my $taxid (@taxidOrder) {
    my %content;
    $content{'id'} = $taxid;
    $content{'metadata'}{'taxonomy'} = $taxidTaxIdx{$taxid}{'tax'};

    push(@newRows,\%content);
  }

  my @newData;
  for (my $iNew=0;$iNew<@taxidOrder;$iNew++) {
    my $taxid = $taxidOrder[$iNew];
    foreach my $i (@{$taxidTaxIdx{$taxid}{'idx'}}) {
      for (my $j=0;$j<@{$jsonHashRef->{'data'}[$i]};$j++) {
        $newData[$iNew][$j] += $jsonHashRef->{'data'}[$i][$j];
      }
    }
  }

  my %outputHash = %{$jsonHashRef};
  $outputHash{'rows'} = \@newRows;
  $outputHash{'data'} = \@newData;
  $outputHash{'shape'}[0] = scalar @newRows;

  &writeBiom(\%outputHash,$output);
}

## task: mergeBiomById
sub mergeBiom { # arg: inputFileNames, outputFileName
  my ($input,$output) = @_;

  if (scalar (@{[split(',',$input)]}) <= 1) {
    print "Invalid number of input file\n";
    &usage();
    exit;
  }
  if (scalar (@{[split(',',$output)]}) > 1) {
    print "Invalid number of input file\n";
    &usage();
    exit;
  }

  my %outputHash;
  $outputHash{'id'} = $output;
  $outputHash{'id'} =~ s/\.biom$//i;

  my @allInput = split(',',$input);

  my %taxInfo;
  my @sampleOrder;
  my %table;
  foreach my $file (@allInput) {
    if (!-e $file) {
      next;
    }
#&checkBiom($file);
#print '--'."\n";

    my $jsonHashRef = &readBiom($file);

    $outputHash{'format'} = $jsonHashRef->{'format'};
    $outputHash{'format_url'} = $jsonHashRef->{'format_url'};
    $outputHash{'type'} = $jsonHashRef->{'type'};
    $outputHash{'generated_by'} = $jsonHashRef->{'generated_by'};
    $outputHash{'date'} = $jsonHashRef->{'date'};
    $outputHash{'matrix_type'} = $jsonHashRef->{'matrix_type'};
    $outputHash{'matrix_element_type'} = $jsonHashRef->{'matrix_element_type'};

    if ($jsonHashRef->{'shape'}[1] != 1) {
      print "Invalid matrix shape: ".encode_json($jsonHashRef->{'shape'})."\n";
      next;
    }

    for (my $i=0;$i<$jsonHashRef->{'shape'}[0];$i++) {
      my $id = $jsonHashRef->{'rows'}[$i]{'id'};
      my $metadata = $jsonHashRef->{'rows'}[$i]{'metadata'};
      my $level = 0;
      for (my $j=0;$j<@{$jsonHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'}};$j++) {
        if ($jsonHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'}[$j] !~ /unclassified/) {
          $level = $j;
        }
      }

      if (!defined($taxInfo{$id}) || ($taxInfo{$id}{'level'} < $level)) {
        $taxInfo{$id}{'metadata'} = $metadata;
        $taxInfo{$id}{'level'} = $level;
      }
    }

    push(@sampleOrder,$jsonHashRef->{'columns'}[0]{'id'});
    push(@{$outputHash{'columns'}},$jsonHashRef->{'columns'}[0]);

    for (my $i=0;$i<$jsonHashRef->{'shape'}[0];$i++) {
      my $taxid = $jsonHashRef->{'rows'}[$i]{'id'};
      my $sample = $jsonHashRef->{'columns'}[0]{'id'};
      my $val = $jsonHashRef->{'data'}[$i][0];

      $table{$taxid}{$sample} += $val;
    }
  }

  my @taxOrder = sort{$a cmp $b} (keys %table);

  if ((scalar @sampleOrder) != (scalar &uniq(@sampleOrder))) {
    print "Duplicated sample ID!\n";
    exit;
  }

  for (my $i=0;$i<@taxOrder;$i++) {
    my %content;
    $content{'id'} = $taxOrder[$i];
    $content{'metadata'} = $taxInfo{$taxOrder[$i]}{'metadata'};

    push(@{$outputHash{'rows'}},\%content);

    for (my $j=0;$j<@sampleOrder;$j++) {
      $outputHash{'data'}[$i][$j] = $table{$taxOrder[$i]}{$sampleOrder[$j]} + 0;
    }
  }

  $outputHash{'shape'} = [(scalar @taxOrder),(scalar @sampleOrder)];

  &writeBiom(\%outputHash,$output);
#&checkBiom($output);
}

## task: splitBiomById
sub splitBiom { # arg: inputFileName
  my $input = shift;
  my $skipZero = shift;

  if (scalar (@{[split(',',$input)]}) > 1) {
    print "Invalid number of input file\n";
    &usage();
    exit;
  }

  my $jsonHashRef = &readBiom($input);

  for (my $j=0;$j<$jsonHashRef->{'shape'}[1];$j++) {
    my $id = $jsonHashRef->{'columns'}[$j]{'id'};

    my %outputHash = %{$jsonHashRef};
    $outputHash{'id'} = $id;

    $outputHash{'columns'} = [$jsonHashRef->{'columns'}[$j]];

    $outputHash{'rows'} = [];
    $outputHash{'data'} = [];
    for (my $i=0;$i<$jsonHashRef->{shape}[0];$i++) {
      my $val = $jsonHashRef->{'data'}[$i][$j];

      if (($skipZero == 1) && ($val == 0)) {
        next;
      }

      push(@{$outputHash{'rows'}},$jsonHashRef->{'rows'}[$i]);
      push(@{$outputHash{'data'}},[$val]);
    }

    $outputHash{'shape'} = [(scalar @{$outputHash{'rows'}}),1];

    &writeBiom(\%outputHash,$id.'_split.biom');
  }
}

## task: sparseToDense
sub sparseToDense { # arg: inputFileName , outputFileName
  my $input = shift;
  my $output = shift;

  $output =~ s/\.biom$//;

  if (scalar (@{[split(',',$input)]}) > 1) {
    print "Invalid number of input file\n";
    &usage();
    exit;
  }

  my $jsonHashRef = &readBiom($input);

#  &writeBiom($jsonHashRef,$output.'.biom');
  if ($jsonHashRef->{'matrix_type'} ne 'sparse') {
    print "Invalid matrix_type in $input\n";
    &usage();
    exit;
  }

  my %outputHash;
  $outputHash{'id'} = $output;
  $outputHash{'format'} = $jsonHashRef->{'format'};
  $outputHash{'format_url'} = $jsonHashRef->{'format_url'};
  $outputHash{'type'} = $jsonHashRef->{'type'};
  $outputHash{'generated_by'} = $jsonHashRef->{'generated_by'};
  $outputHash{'date'} = $jsonHashRef->{'date'};
  $outputHash{'matrix_type'} = 'dense';
  $outputHash{'matrix_element_type'} = 'int';

  $outputHash{'rows'} = $jsonHashRef->{'rows'};
  $outputHash{'columns'} = $jsonHashRef->{'columns'};
  $outputHash{'shape'} = $jsonHashRef->{'shape'};

  for (my $i=0;$i<$jsonHashRef->{'shape'}[0];$i++) {
    for (my $j=0;$j<$jsonHashRef->{'shape'}[1];$j++) {
      $outputHash{'data'}[$i][$j] = 0;
    }
  }

  foreach my $row (@{$jsonHashRef->{'data'}}) {
    my ($i,$j,$val) = @{$row};

    $outputHash{'data'}[$i][$j] = $val;
  }

  &writeBiom(\%outputHash,$output.'.biom');
#"columns"
#"data"
#"shape"
#"rows"

=h
  for (my $j=0;$j<$jsonHashRef->{'shape'}[1];$j++) {
    my $id = $jsonHashRef->{'columns'}[$j]{'id'};

    my %outputHash = %{$jsonHashRef};
    $outputHash{'id'} = $id;

    $outputHash{'columns'} = [$jsonHashRef->{'columns'}[$j]];

    $outputHash{'rows'} = [];
    $outputHash{'data'} = [];
    for (my $i=0;$i<$jsonHashRef->{shape}[0];$i++) {
      my $val = $jsonHashRef->{'data'}[$i][$j];

      if (($skipZero == 1) && ($val == 0)) {
        next;
      }

      push(@{$outputHash{'rows'}},[$jsonHashRef->{'rows'}[$i]]);
      push(@{$outputHash{'data'}},[$val]);
    }

    $outputHash{'shape'} = [(scalar @{$outputHash{'rows'}}),1];

    &writeBiom(\%outputHash,$id.'_split.biom');
  }
=cut
}

## task: idChange
sub idChange { # arg: inputFileName, outputFileName
  my ($input,$output,$names) = @_;

  if (scalar (@{[split(',',$input)]}) > 1) {
    print "Invalid number of input file\n";
    &usage();
    exit;
  }
  if (scalar (@{[split(',',$output)]}) > 1) {
    print "Invalid number of input file\n";
    &usage();
    exit;
  }
  if (scalar (@{[split(',',$names)]}) > 1) {
    print "Invalid number of input file\n";
    &usage();
    exit;
  }

  my $jsonHashRef = &readBiom($input);

  if (!-e $names) {
    print "File $names not found!\n";
    &usage();
    exit;
  }

  my %idMapping;
  open(IN,'<'.$names);
  while (<IN>) {
    chomp;
    if (($_ eq '') || ($_ =~ /^#/)) {
      next;
    }

    my ($from,$to) = split("\t",$_);

    if (defined($idMapping{$from}) && ($idMapping{$from} ne '') && ($idMapping{$from} ne $to)) {
      print "Multiple result with id:$from\n";
      exit;
    }

    $idMapping{$from} = $to;
  }
  close(IN);

  for (my $i=0;$i<$jsonHashRef->{'shape'}[0];$i++) {
    if (!defined($jsonHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'})) {
      next;
    }

    my $id = $jsonHashRef->{'rows'}[$i]{'id'};

    if (defined($idMapping{$id}) && ($idMapping{$id} ne '')) {
      $jsonHashRef->{'rows'}[$i]{'id'} = $idMapping{$id};
    }
  }

  &writeBiom($jsonHashRef,$output);
}

sub usage {
  print <<USAGE;
Usage: perl $0 task inputFileName [outputFileName] [idMappingFile]

 - uniqByTax
    Combine rows with same taxonomy information from biom file.
    Required parameter:
      inputFileName: single filename
      outputFileName: single filename
 - uniqById
    Combine rows with same ID (TaxID) from biom file.
    Required parameter:
      inputFileName: single filename
      outputFileName: single filename
 - mergeBiom
    Merge multiple biom files to one biom file.
    Required parameter:
      inputFileName: comma-separate multiple biom filename
      outputFileName: single filename
 - splitBiom
    Split one biom file to multiple biom, use column ID (sample name) as filename.
    Required parameter:
      inputFileName: single filename
 - splitBiomSkip0
    Split one biom file to multiple biom, use column ID (sample name) as filename. Filter rows (tax) with zero count.
    Required parameter:
      inputFileName: single filename
 - sparseToDense
    Format modify from sparse to dense
    Required parameter:
      inputFileName: single biom filename in sparse
      outputFileName: single biom filename
 - idChange
    Modify row id using idMappingFile
    Required parameter:
      inputFileName: single biom filename
      outputFileName: single biom filename
      idMappingFile: id mapping table with tab-separated two column id
 - idChangeUsingFasta
    Modify row id using idMappingFile
    Required parameter:
      inputFileName: single biom filename
      outputFileName: single biom filename
      idMappingFile: two fileName separated by comma. 1. .rep.fasta  2. id mapping table with tab-separated two column id
USAGE
}

if (@ARGV < 2) {
  &usage();
  exit;
}

my $task = shift;
my $input = shift;
my $output = shift;
my $file = shift;

if (!defined($task) || !defined($input)) {
  &usage();
  exit;
}

if (lc($task) eq 'uniqbytax') {
  if (!defined($output) || ($output eq '')) {
    &usage();
      exit;
  }

  &uniqByTax($input,$output);
}
elsif (lc($task) eq 'uniqbyid') {
  if (!defined($output) || ($output eq '')) {
    &usage();
      exit;
  }

  &uniqById($input,$output);
}
elsif (lc($task) eq 'mergebiom') {
  if (!defined($output) || ($output eq '')) {
    &usage();
      exit;
  }

  &mergeBiom($input,$output);
}
elsif (lc(substr($task,0,9)) eq 'splitbiom') {
  my $skip = 0;
  if ($task =~ /skip0/i) {
    $skip = 1;
  }
  &splitBiom($input,$skip);
}
elsif (lc($task) eq 'idchange') {
  if (!defined($file) || ($file eq '')) {
    &usage();
      exit;
  }

  &idChange($input,$output,$file);
}
elsif (lc($task) eq 'idchangeusingfasta') {
  my ($fasta,$tab) = split(',',$file);

  my %fastaName;
  open(IN,'grep "^>" '.$fasta.' | ');
  while (<IN>) {
    chomp;
    if ($_ eq '') {
      next;
    }

    my $name = substr($_,1);
    ($name) = split(/\|/,$name);
    my ($to,$from) = split(/\s+/,$name);

    $fastaName{$from} = $to;
  }
  close(IN);

  my %tabName;
  open(IN,'<'.$tab);
  while (<IN>) {
    chomp;
    if (($_ eq '') || ($_ =~ /^#/)) {
      next;
    }

    my ($from,$to) = split("\t",$_);

    $tabName{$from} = $to;
  }
  close(IN);

  open(OUT,'>tmp.tab');
  foreach my $aa (keys %fastaName) {
    if (!defined($fastaName{$aa}) || ($fastaName{$aa} eq '')) {
      next;
    }

    my $bb = $fastaName{$aa};
    if (!defined($tabName{$bb}) || ($tabName{$bb} eq '')) {
      next;
    }

    my $cc = $tabName{$bb};

    print OUT join("\t",$aa,$cc)."\n";
  }
  close(OUT);

  &idChange($input,$output,'tmp.tab');

  unlink('tmp.tab');
}
elsif (lc($task) eq 'sparsetodense') {
  &sparseToDense($input,$output);
}
elsif (lc($task) eq 'checkbiom') {
  &checkBiom($input);
}
else {
  &usage();
  exit
}
