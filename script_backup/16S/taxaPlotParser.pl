#!/usr/bin/perl
use strict;

if ((@ARGV == 0) || !(-e $ARGV[0].'/otu_table_L3.txt')) {
  print "perl $0 [taxaplot raw data] [output folder]\n";
  exit;
}

my $input = shift(@ARGV);
my $output = shift(@ARGV);

$input =~ s/\/+$//;
$output =~ s/\/+$//;

if (!$output) {
  $output = '.';
}

use Spreadsheet::WriteExcel;
my $workbook = Spreadsheet::WriteExcel->new($output.'/taxaPlot_data.xls');

my $centerAlignHeader1 = $workbook->add_format(valign=>'vcenter',align=>'center',bold=>1);
my $centerAlignHeader2 = $workbook->add_format(align=>'center',bold=>1);
my $centerAlign = $workbook->add_format(valign=>'vcenter',align=>'center');

my @level = qw / Kindom Phylum Class Order Family Genus Species / ;
unshift(@level,'');
unshift(@level,'');

foreach my $id (2..8) {
  if (!-e $input.'/otu_table_L'.$id.'.txt') {
    next;
  }

  my @taxOrder;

  my %readTable;
  open(IN,'<'.$input.'/otu_table_L'.$id.'.txt');
  chomp(my $header = <IN>);
  my @colNames = split("\t",$header);
  while (<IN>) {
    chomp;
    if ($_ eq '') {
      next;
    }

    my @array = split("\t",$_);

    my $rowName = $array[0];
    push(@taxOrder,$rowName);

    for (my $i=1;$i<@array;$i++) {
      my $colName = $colNames[$i];
      my $val = $array[$i];

      $readTable{$rowName}{$colName} += $val;
      $readTable{'Sum'}{$colName} += $val;
      $readTable{$rowName}{'Sum'} += $val;
      $readTable{'Sum'}{'Sum'} += $val;
    }
  }
  close(IN);

  my $worksheet = $workbook->add_worksheet($level[$id]);

  $worksheet->merge_range(0,0,1,0,'#OTU ID',$centerAlignHeader1);
  $worksheet->merge_range(0,1,0,2,'Total',$centerAlignHeader1);
  $worksheet->write(1,1,'Count',$centerAlignHeader2);
  $worksheet->write(1,2,'%',$centerAlignHeader2);

  for (my $i=1;$i<@colNames;$i++) {
    $worksheet->merge_range(0,(2*$i+1),0,(2*$i+2),$colNames[$i],$centerAlignHeader1);
    $worksheet->write(1,(2*$i+1),'Count',$centerAlignHeader2);
    $worksheet->write(1,(2*$i+2),'%',$centerAlignHeader2);
  }

  my $row = 2;
  foreach my $rowName (@taxOrder) {
    $worksheet->write($row,0,$rowName);
    $worksheet->write_number($row,1,$readTable{$rowName}{'Sum'});
    $worksheet->write_number($row,2,(100*$readTable{$rowName}{'Sum'}/$readTable{'Sum'}{'Sum'}));

    for (my $i=1;$i<@colNames;$i++) {
      $worksheet->write_number($row,(2*$i+1),$readTable{$rowName}{$colNames[$i]});
      $worksheet->write_number($row,(2*$i+2),(100*$readTable{$rowName}{$colNames[$i]}/$readTable{'Sum'}{$colNames[$i]}));
    }

    $row ++;
  }

  $worksheet->freeze_panes(2,1);
}

$workbook->close();










=h
otu_table_L3.txt
Taxon	C1	T0	T3.1.2	T3.1	T4	T5a	T5b	T5c
Root;k__Archaea;p__Crenarchaeota	55.0	160.0	302.0	249.0	52.0	174.0	286.0	124.0
Root;k__Archaea;p__Euryarchaeota	26.0	90.0	87.0	66.0	21.0	303.0	517.0	45.0
Root;k__Archaea;p__[Parvarchaeota]	10.0	38.0	61.0	47.0	10.0	37.0	23.0	22.0
Root;k__Archaea;unclassified	1.0	5.0	4.0	6.0	2.0	3.0	12.0	4.0
Root;k__Bacteria;p__AC1	57.0	70.0	52.0	64.0	42.0	48.0	77.0	83.0
Root;k__Bacteria;p__AD3	0.0	1.0	0.0	0.0	1.0	0.0	0.0	0.0
Root;k__Bacteria;p__Acidobacteria	7676.0	7799.0	6967.0	6141.0	4830.0	6469.0	10180.0	8329.0
Root;k__Bacteria;p__Actinobacteria	4062.0	4200.0	3684.0	3290.0	1769.0	3347.0	3837.0	3373.0
Root;k__Bacteria;p__Aquificae	1.0	0.0	0.0	1.0	0.0	0.0	0.0	0.0
=cut
