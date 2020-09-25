#!/usr/bin/env perl
use strict;

sub uniq {
  my %seen;
  map $seen{$_}=1,@_;
  my @uniq = sort{$a <=> $b} (keys %seen);
  return(@uniq);
}

#my @matrix = qw / ACE chao1 enspie fisher_alpha goods_coverage observed_species shannon simpson_reciprocal / ;
my @matrix = qw / chao1 enspie fisher_alpha goods_coverage observed_species shannon simpson_reciprocal / ;

my $input = shift(@ARGV);
my $output = shift(@ARGV);

$input =~ s/\/+$//g;
$output =~ s/\/+$//g;

if (!$output) {
  $output = '.';
}

{
if (($input eq '') || (!-e $input) || (!-d $input)) {
  print "Usage: perl $0 [rare step 3 folder] [output folder]\n";
  exit;
}
my $check = 1;
foreach my $mat (@matrix) {
  if (!-e $input.'/'.$mat.'.txt') {
    $check = 0;
    last;
  }
}
if ($check == 0) {
  print "Something in rare step 3 folder missing\n";
  exit;
}
}

my @allSPS;
my @sampleOrder;
my %readTable;
foreach my $mat (@matrix) {
  my $file = $input.'/'.$mat.'.txt';
  open(IN,'<'.$file);
  chomp(my $header = <IN>);
  my @samples = split("\t",$header);
  if (!defined($sampleOrder[0])) {
    @sampleOrder = @samples;
  }
  while (<IN>) {
    chomp;
    if ($_ eq '') {
      next;
    }

    my @array = split("\t",$_);

    my $seqPerSample = $array[1];
    push(@allSPS,$seqPerSample);
    for (my $i=3;$i<@array;$i++) {
      my $sampleName = $samples[$i];
      my $val = $array[$i];

      if (($val eq 'n/a') || ($val eq 'nan')) {
        next;
      }

      #$readTable{$mat}{$sampleName}{$seqPerSample} = $val;
      push(@{$readTable{$mat}{$sampleName}{$seqPerSample}{'values'}},$val);
      $readTable{$mat}{$sampleName}{$seqPerSample}{'sum'} += $val;
      $readTable{$mat}{$sampleName}{$seqPerSample}{'n'} ++;
    }
  }
  close(IN);
}

shift(@sampleOrder);
shift(@sampleOrder);
shift(@sampleOrder);

@allSPS = &uniq(@allSPS);

use Spreadsheet::WriteExcel;
my $workbook = Spreadsheet::WriteExcel->new($output.'/rarePlot_data.xls');
my $centerAlignHeader1 = $workbook->add_format(valign=>'vcenter',align=>'center',bold=>1);
my $centerAlignHeader2 = $workbook->add_format(align=>'center',bold=>1);
my $centerAlign = $workbook->add_format(valign=>'vcenter',align=>'center');

foreach my $mat (@matrix) {
  my $worksheet = $workbook->add_worksheet($mat);

  $worksheet->merge_range(0,0,1,0,'Sample name',$centerAlignHeader1);
  $worksheet->merge_range(0,1,0,($#allSPS+1),'Rarefaction measurement on different number of sequences per sample',$centerAlignHeader1);

  for (my $i=0;$i<@allSPS;$i++) {
    $worksheet->write(1,($i+1),$allSPS[$i],$centerAlignHeader2);
  }

  my $row = 2;
  foreach my $sampleName (@sampleOrder) {
    $worksheet->write($row,0,$sampleName);

    for (my $i=0;$i<@allSPS;$i++) {
      if (defined($readTable{$mat}{$sampleName}{$allSPS[$i]}) && (lc($readTable{$mat}{$sampleName}{$allSPS[$i]}) ne 'nan')) {
        my $val = $readTable{$mat}{$sampleName}{$allSPS[$i]}{'sum'} / $readTable{$mat}{$sampleName}{$allSPS[$i]}{'n'};
        $val += 0;

        #$worksheet->write_number($row,($i+1),$readTable{$mat}{$sampleName}{$allSPS[$i]});
        $worksheet->write_number($row,($i+1),$val);
      }
      else {
        $worksheet->write($row,($i+1),'NaN');
      }
    }

    $row ++;
  }

  $worksheet->freeze_panes(2,1);
}

$workbook->close();
