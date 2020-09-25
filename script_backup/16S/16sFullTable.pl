#!/usr/bin/env perl
use strict;

if ((!-e 'B.5_OTU_cluster/') || (!-e 'B.5_OTU_cluster/97/')) {
  exit;
}

my @samples;

my %allInfo;
my @header;
my %totalCount;
my %mergeTax;
open(IN,'<B.5_OTU_cluster/otu_table.taxID.filter.dense.dedup.txt');
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  my @array = split("\t",$_);

  if ($_ =~ /^#/) {
    if (@array < 3) {
      next;
    }

    @header = @array;

    @samples = @header[1..($#array-1)];

    next;
  }

  my $otuID = $array[0];
  my $tax = $array[-1];

  for (my $i=1;$i<$#array;$i++) {
    $allInfo{$otuID}{'count'}{$header[$i]} = $array[$i];
    $totalCount{$header[$i]} += $array[$i];
  }

  $mergeTax{$otuID} = $tax;
}
close(IN);

foreach my $sn (@samples) {
  my $file = 'B.5_OTU_cluster/97/'.$sn.'.otu_table.taxID.filter.dense.dedup.txt';
  if (!-e $file) {
    next;
  }

  open(IN,'<'.$file);
  while (<IN>) {
    chomp;
    if (($_ eq '') || ($_ =~ /^#/)) {
      next;
    }

    my @array = split("\t",$_);

    if (!defined($allInfo{$array[0]}{'count'}{$sn}) || ($allInfo{$array[0]}{'count'}{$sn} != $array[1])) {
      print 'Unequal: '.$array[0]."\n";
      next;
    }

    $allInfo{$array[0]}{'taxonomy'}{$sn} = $array[2];
  }
  close(IN);
}

if (!-e 'report_html/data') {
  system('mkdir -p report_html/data/');
}

use Excel::Writer::XLSX;
my $workbook = Excel::Writer::XLSX->new('report_html/data/FullTable.otu.xlsx');

my $centerAlignHeader = $workbook->add_format(align=>'center',bold=>1);

foreach my $i ('count','taxonomy','ratio') {
  my $worksheet = $workbook->add_worksheet( ($i eq 'ratio')? $i.' (%)' : $i );
  $worksheet->write(0,0,'OTU ID',$centerAlignHeader);
  $worksheet->write(0,1,'Taxonomy',$centerAlignHeader);
  for (my $c=0;$c<@samples;$c++) {
    $worksheet->write_string(0,$c+2,$samples[$c],$centerAlignHeader);
  }

  my $nl = 1;
  open(OUT,'>FullTable.otu.'.$i.'.txt');
  print OUT join("\t",'OTU ID','Taxonomy',@samples)."\n";
  foreach my $otuID (sort{$a cmp $b} (keys %allInfo)) {
    $worksheet->write_string($nl,0,$otuID);
    $worksheet->write_string($nl,1,$mergeTax{$otuID});

    print OUT $otuID."\t".$mergeTax{$otuID};

    for (my $c=0;$c<@samples;$c++) {
      my $sn = $samples[$c];

      my $val = ($i eq 'taxonomy') ? 'NA' : 0;
      if (defined($allInfo{$otuID}{$i}{$sn})) {
        $val = $allInfo{$otuID}{$i}{$sn};
      }
      if ($i eq 'ratio') {
        $val = $allInfo{$otuID}{'count'}{$sn} + 0;
        $val = $val*100 / $totalCount{$sn};
      }
      if ($i eq 'taxonomy') {
        if ($val eq $mergeTax{$otuID}) {
          $val = '--';
        }
      }

      if ($i eq 'taxonomy') {
        $worksheet->write_string($nl,$c+2,$val);
      }
      else {
        $worksheet->write_number($nl,$c+2,$val);
      }
      print OUT "\t".$val;
    }

    $nl ++;
    print OUT "\n";

    $worksheet->freeze_panes(1,2);
  }
  close(OUT);
}

$workbook->close();

