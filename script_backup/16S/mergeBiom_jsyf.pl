#!/usr/bin/perl
use strict;

use JSON;

=h
my $jsonStr = '';
while (<>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  $jsonStr .= $_;
}


my $biomStructure  = decode_json $jsonStr;

#print join("\t",(sort{$a cmp $b} (keys %{$biomStructure})))."\n";

my $count = scalar @{$biomStructure->{'rows'}};

for (my $i=0;$i<$count;$i++) {
  my $id = $biomStructure->{'rows'}[$i]{'id'};
  my $val = $biomStructure->{'data'}[$i][0];

  print join("\t",$id,$val)."\n";
}

=cut

#combine biom files

my $outputName = shift(@ARGV);
if ((!$outputName) || ($outputName =~ /\//)) { $outputName = 'Merged';}

my @files;
foreach my $file (@ARGV) {
  if (($file ne '') && (-e $file)) {
    push(@files,$file);
  }
}

my %combineBiom;
$combineBiom{'id'} = $outputName;
$combineBiom{'format'} = 'Biological Observation Matrix 0.9.1';
$combineBiom{'format_url'} = 'http://biom-format.org';
$combineBiom{'type'} = 'OTU table';
$combineBiom{'generated_by'} = 'mothur1.33.3';
$combineBiom{'date'} = 'Fri Jun 30 05:58:38 2017';
$combineBiom{'columns'} = [{"id" => $outputName, "metadata" => undef}];
$combineBiom{'matrix_type'} = 'dense';
$combineBiom{'matrix_element_type'} = 'int';

my %biomRows;
foreach my $file (@files) {
  open(IN,'<'.$file);
  chomp(my @content = <IN>);
  my $thisBiom = decode_json(join('',@content));
  close(IN);

  $combineBiom{'date'} = $thisBiom->{'date'};

  my $n = scalar @{$thisBiom->{'rows'}};

  for (my $i=0;$i<$n;$i++) {
    my $id = $thisBiom->{'rows'}[$i]{'id'};
    my $metadata = $thisBiom->{'rows'}[$i]{'metadata'};
    my $data = $thisBiom->{'data'}[$i][0];

    $biomRows{$id}{'metadata'} = $metadata;
    $biomRows{$id}{'data'} += $data;
  }
}

my $nRows = scalar (keys %biomRows);

foreach my $id (sort{$a cmp $b} (keys %biomRows)) {
  push(@{$combineBiom{'rows'}},{'id' => "$id",'metadata' => $biomRows{$id}{'metadata'}});
  push(@{$combineBiom{'data'}},[$biomRows{$id}{'data'}]);
}

$combineBiom{'shape'} = [$nRows,1];

open(OUT,'>'.$outputName.'.biom');
my $output = encode_json(\%combineBiom);
print OUT $output."\n";
close(OUT);
