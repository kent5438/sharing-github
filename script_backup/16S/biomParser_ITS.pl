#!/usr/bin/perl
use strict;

if (@ARGV < 2) {
  print "Usage: perl $0 inputBiomFile outputBiomFile\n";
  exit;
}

use JSON;
my $biomJson = '';
open(IN,'<'.$ARGV[0]);
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  $biomJson .= $_;
}
close(IN);

my $biomHashRef = decode_json($biomJson);
my %combineRows;
for (my $i=0;$i<@{$biomHashRef->{'rows'}};$i++) {
  my $otuID = $biomHashRef->{'rows'}[$i]{'id'};
  my @taxonomy;
  if (ref($biomHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'}) eq 'ARRAY') {
    @taxonomy = @{$biomHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'}};
  }
  else {
    @taxonomy = split('; ',$biomHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'});
  }

  if ($taxonomy[0] !~ /Root/) {
    unshift(@taxonomy,'Root');
  }

  my @bootstrap;
  if (defined($biomHashRef->{'rows'}[$i]{'metadata'}{'bootstrap'})) {
    @bootstrap = @{$biomHashRef->{'rows'}[$i]{'metadata'}{'bootstrap'}};
  }

  while ($#bootstrap < $#taxonomy) {
    unshift(@bootstrap,100);
  }

  $combineRows{$taxonomy[0].':'.$bootstrap[0]}{$taxonomy[1].':'.$bootstrap[1]}{$taxonomy[2].':'.$bootstrap[2]}{$taxonomy[3].':'.$bootstrap[3]}{$taxonomy[4].':'.$bootstrap[4]}{$taxonomy[5].':'.$bootstrap[5]}{$taxonomy[6].':'.$bootstrap[6]}{$taxonomy[7].':'.$bootstrap[7]}{'taxonomy'} = \@taxonomy;
  $combineRows{$taxonomy[0].':'.$bootstrap[0]}{$taxonomy[1].':'.$bootstrap[1]}{$taxonomy[2].':'.$bootstrap[2]}{$taxonomy[3].':'.$bootstrap[3]}{$taxonomy[4].':'.$bootstrap[4]}{$taxonomy[5].':'.$bootstrap[5]}{$taxonomy[6].':'.$bootstrap[6]}{$taxonomy[7].':'.$bootstrap[7]}{'bootstrap'} = \@bootstrap;
  push(@{$combineRows{$taxonomy[0].':'.$bootstrap[0]}{$taxonomy[1].':'.$bootstrap[1]}{$taxonomy[2].':'.$bootstrap[2]}{$taxonomy[3].':'.$bootstrap[3]}{$taxonomy[4].':'.$bootstrap[4]}{$taxonomy[5].':'.$bootstrap[5]}{$taxonomy[6].':'.$bootstrap[6]}{$taxonomy[7].':'.$bootstrap[7]}{'otuID'}},$otuID);
  $combineRows{$taxonomy[0].':'.$bootstrap[0]}{$taxonomy[1].':'.$bootstrap[1]}{$taxonomy[2].':'.$bootstrap[2]}{$taxonomy[3].':'.$bootstrap[3]}{$taxonomy[4].':'.$bootstrap[4]}{$taxonomy[5].':'.$bootstrap[5]}{$taxonomy[6].':'.$bootstrap[6]}{$taxonomy[7].':'.$bootstrap[7]}{'count'} += $biomHashRef->{'data'}[$i][0];
}

$biomHashRef->{'rows'} = [];
$biomHashRef->{'data'} = [];

foreach my $root (sort{$a cmp $b} (keys %combineRows)) {
  foreach my $kk (sort{$a cmp $b} (keys %{$combineRows{$root}})) {
    foreach my $pp (sort{$a cmp $b} (keys %{$combineRows{$root}{$kk}})) {
      foreach my $cc (sort{$a cmp $b} (keys %{$combineRows{$root}{$kk}{$pp}})) {
        foreach my $oo (sort{$a cmp $b} (keys %{$combineRows{$root}{$kk}{$pp}{$cc}})) {
          foreach my $ff (sort{$a cmp $b} (keys %{$combineRows{$root}{$kk}{$pp}{$cc}{$oo}})) {
            foreach my $gg (sort{$a cmp $b} (keys %{$combineRows{$root}{$kk}{$pp}{$cc}{$oo}{$ff}})) {
              foreach my $ss (sort{$a cmp $b} (keys %{$combineRows{$root}{$kk}{$pp}{$cc}{$oo}{$ff}{$gg}})) {
                my %outputRow;
                $outputRow{'id'} = $combineRows{$root}{$kk}{$pp}{$cc}{$oo}{$ff}{$gg}{$ss}{'otuID'}[0];
                $outputRow{'metadata'}{'taxonomy'} = $combineRows{$root}{$kk}{$pp}{$cc}{$oo}{$ff}{$gg}{$ss}{'taxonomy'};
                $outputRow{'metadata'}{'bootstrap'} = $combineRows{$root}{$kk}{$pp}{$cc}{$oo}{$ff}{$gg}{$ss}{'bootstrap'};

                push(@{$biomHashRef->{'rows'}},\%outputRow);

                push(@{$biomHashRef->{'data'}},[$combineRows{$root}{$kk}{$pp}{$cc}{$oo}{$ff}{$gg}{$ss}{'count'}]);
              }
            }
          }
        }
      }
    }
  }
}

my $nRow = scalar @{$biomHashRef->{'rows'}};
my $nData = scalar @{$biomHashRef->{'data'}};

if ($nRow != $nData) {
  print "Uneqal data length\n";
  exit;
}

$biomHashRef->{'shape'}[0] = $nData;

open(OUT,'>'.$ARGV[1]);
#my $output = encode_json($biomHashRef);
#print OUT $output."\n";
print OUT '{'."\n";
print OUT '      "id":"'.$biomHashRef->{'id'}.'",'."\n";
print OUT '      "format": "'.$biomHashRef->{'format'}.'",'."\n";
print OUT '      "format_url": "'.$biomHashRef->{'format_url'}.'",'."\n";
print OUT '      "type": "'.$biomHashRef->{'type'}.'",'."\n";
print OUT '      "generated_by": "'.$biomHashRef->{'generated_by'}.'",'."\n";
print OUT '      "date": "'.$biomHashRef->{'date'}.'",'."\n";
print OUT '      "rows":['."\n";
for (my $i=0;$i<$nRow;$i++) {
  if ($i == ($nRow-1)) {
    print OUT '            {"id":"'.$biomHashRef->{'rows'}[$i]{'id'}.'", "metadata":{"taxonomy":["'.join('", "',@{$biomHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'}}).'"], "bootstrap":['.join(', ',@{$biomHashRef->{'rows'}[$i]{'metadata'}{'bootstrap'}}).']}}'."\n";
  }
  else {
    print OUT '            {"id":"'.$biomHashRef->{'rows'}[$i]{'id'}.'", "metadata":{"taxonomy":["'.join('", "',@{$biomHashRef->{'rows'}[$i]{'metadata'}{'taxonomy'}}).'"], "bootstrap":['.join(', ',@{$biomHashRef->{'rows'}[$i]{'metadata'}{'bootstrap'}}).']}},'."\n";
  }
}
print OUT '      ],'."\n";
print OUT '      "columns":['."\n";
print OUT '            {"id":"'.$biomHashRef->{'columns'}[0]{'id'}.'", "metadata":'.(defined($biomHashRef->{'columns'}[0]{'metadata'})?'"'.$biomHashRef->{'columns'}[0]{'metadata'}.'"':'null').'}'."\n";
print OUT '      ],'."\n";
print OUT '      "matrix_type": "'.$biomHashRef->{'matrix_type'}.'",'."\n";
print OUT '      "matrix_element_type": "'.$biomHashRef->{'matrix_element_type'}.'",'."\n";
print OUT '      "shape": ['.join(',',@{$biomHashRef->{'shape'}}).'],'."\n";
for (my $i=0;$i<$nData;$i++) {
  if ($i == 0) {
    print OUT '      "data":  [['.$biomHashRef->{'data'}[$i][0].'],'."\n";
  }
  elsif ($i == ($nData-1)) {
    print OUT '            ['.$biomHashRef->{'data'}[$i][0].']]'."\n";
  }
  else {
    print OUT '            ['.$biomHashRef->{'data'}[$i][0].'],'."\n";
  }
}
print OUT '}'."\n";
close(OUT);
