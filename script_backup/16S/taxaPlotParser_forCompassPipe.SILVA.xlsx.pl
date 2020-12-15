#!/usr/bin/perl
use strict;

if ((@ARGV == 0) || !(-e $ARGV[0].'/otu_table.taxID.filter.dense.dedup_L3.txt')) {
  print "perl $0 [taxaplot raw data] [output folder]\n";
  exit;
}

my $input = shift(@ARGV);
my $output = (@ARGV > 0)?shift(@ARGV):'.';

$input .= '/';
$output .= '/';

$input =~ s/\/+/\//;
$output =~ s/\/+/\//;

my @level = qw / Kingdom Phylum Class Order Family Genus Species / ;
#unshift(@level,'');
unshift(@level,'');

my %levelTable;
foreach my $levelId (1..7) {
  my $file = $input.'otu_table.taxID.filter.dense.dedup_L'.$levelId.'.txt';
  if (!-e $file) {
    next;
  }

  my @colOrder;
  open(IN,'<'.$file);
  while (<IN>) {
    chomp;
    if ($_ eq '') {
      next;
    }

    my @array = split("\t",$_);

    my $tax = shift(@array);

    if ($. == 1) {
      $levelTable{$level[$levelId]}{'columns'} = \@array;
      @colOrder = @array;
      next;
    }

    push(@{$levelTable{$level[$levelId]}{'rows'}},$tax);

    for (my $i=0;$i<@array;$i++) {
      my $val = $array[$i];

      $levelTable{$level[$levelId]}{'data'}{$tax}{$colOrder[$i]} += $val;
      $levelTable{$level[$levelId]}{'data'}{'colSum'}{$colOrder[$i]} += $val;
      $levelTable{$level[$levelId]}{'data'}{$tax}{'rowSum'} += $val;
      $levelTable{$level[$levelId]}{'data'}{'colSum'}{'rowSum'} += $val;
    }
  }
  close(IN);
}

use Excel::Writer::XLSX;
my $workbook = Excel::Writer::XLSX->new($output.'taxaPlot_data.xlsx');

my $centerAlignHeader1 = $workbook->add_format(valign=>'vcenter',align=>'center',bold=>1);
my $centerAlignHeader2 = $workbook->add_format(align=>'center',bold=>1);
my $centerAlign = $workbook->add_format(valign=>'vcenter',align=>'center');

foreach my $dataType ('count','%') {
  foreach my $l (@level) {
    if ($l eq '') {
      next;
    }

    my $worksheet = $workbook->add_worksheet( ($dataType eq 'count')? $l : $l.' ('.$dataType.')' );
    $worksheet->write(0,0,'Taxon',$centerAlignHeader2);
    $worksheet->write(0,1,'Total',$centerAlignHeader2);
    for (my $col=-1;$col<@{$levelTable{$l}{'columns'}};$col++) {
      my $colIdx = $col + 2;

      for (my $row=-1;$row<@{$levelTable{$l}{'rows'}};$row++) {
        my $rowIdx = $row + 1;

        if (($row == -1) && ($col == -1)) {
          next;
        }
        elsif ($row == -1) {
          $worksheet->write(0,$colIdx,$levelTable{$l}{'columns'}[$col],$centerAlignHeader2);
        }
        elsif ($col == -1) {
          $worksheet->write($rowIdx,0,$levelTable{$l}{'rows'}[$row]);

          my $val = $levelTable{$l}{'data'}{$levelTable{$l}{'rows'}[$row]}{'rowSum'} + 0;
          if ($dataType ne 'count') {
            $val = (100 * $val) / $levelTable{$l}{'data'}{'colSum'}{'rowSum'};
          }
          $worksheet->write_number($rowIdx,1,$val);
        }
        else {
          my $colName = $levelTable{$l}{'columns'}[$col];
          my $rowName = $levelTable{$l}{'rows'}[$row];
          my $val = $levelTable{$l}{'data'}{$rowName}{$colName} + 0;

          if ($dataType ne 'count') {
            $val = (100 * $val) / $levelTable{$l}{'data'}{'colSum'}{$colName};
          }

          $worksheet->write_number($rowIdx,$colIdx,$val);
        }
      }
    }

    $worksheet->freeze_panes(1,1);
  }
}

$workbook->close();










=h
otu_table.taxID.filter.dense.dedup_L3.txt
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
