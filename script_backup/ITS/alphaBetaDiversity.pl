#!/usr/bin/perl
use strict;

my $qiime_docker = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` bwawrik/qiime:latest";

sub usageAndExit {
  print "Usage: perl $0 input.biom outputFolder tree\n";
  exit;
}

if (@ARGV < 2) {
  &usageAndExit();
}

my ($input,$output,$tree) = @ARGV;
my $pwd = `pwd`; chomp($pwd);

if($pwd =~ /ITS/){
	if (!$tree) {
		$tree = '/export/EC1680U/metagenomics/ITSdb/align/UNITEv6_sh_99.aln.tre';
	}
	system("rsync -avq $tree .");
}
else { #if($pwd =~ /16S/)
	if (!$tree) {
		$tree = '/export/EC1680U/metagenomics/16Sdb/greenGenes/13_5/gg_13_8_otus/trees/99_otus.tree';
	}
	system("rsync -avq $tree .");
}

$output .= '/';
$output =~ s/\/+/\//g;

if (!(-e $input) || !(-f $input) || !(-e $output) || !(-d $output) || !(-e $tree) || !(-f $tree)) {
  &usageAndExit();
}

#get otuput file name
my $outputFileName = $input;
$outputFileName = pop(@{[split('/',$outputFileName)]});

$outputFileName =~ s/\.\w+$//;
$outputFileName .= '.txt';

# Alpha diversity
if (fork() == 0) {
  my $cmd = $qiime_docker.' alpha_diversity.py -i '.$input.' -o '.$output.'alpha_'.$outputFileName.' -m chao1,enspie,fisher_alpha,goods_coverage,observed_species,shannon,simpson_reciprocal -t 99_otus.tree';
  if ($pwd =~ /ITS/) {
    $cmd = $qiime_docker.' alpha_diversity.py -i '.$input.' -o '.$output.'alpha_'.$outputFileName.' -m chao1,enspie,fisher_alpha,goods_coverage,observed_species,shannon,simpson_reciprocal -t UNITEv6_sh_99.aln.tre';
  }
  print $cmd."\n";
  system($cmd);

  exit;
}

# Beta diversity
if (fork() == 0) {
  my $cmd = $qiime_docker .' beta_diversity.py -i '.$input.' -o '.$output.' -m weighted_unifrac,unweighted_unifrac -t 99_otus.tree';
  if ($pwd =~ /ITS/) {
    $cmd = $qiime_docker .' beta_diversity.py -i '.$input.' -o '.$output.' -m euclidean,pearson -t UNITEv6_sh_99.aln.tre';
  }
  print $cmd."\n";
  system($cmd);

  exit;
}

wait for 0 .. 2;

# Output .xls file of alpha diversity
use Spreadsheet::WriteExcel;
if (fork() == 0) {
  my $workbook = Spreadsheet::WriteExcel->new($output.'alphaDiversity.xls');
  my $centerAlignHeader1 = $workbook->add_format(valign=>'vcenter',align=>'center',bold=>1);
  my $centerAlignHeader2 = $workbook->add_format(align=>'center',bold=>1);
  my $centerAlign = $workbook->add_format(valign=>'vcenter',align=>'center');

  my $worksheet = $workbook->add_worksheet('alpha');
  my $nl = 0;
  open(IN,'<'.$output.'alpha_'.$outputFileName);
  while (<IN>) {
    chomp;
    if ($_ eq '') {
      next;
    }

    my @array = split("\t",$_);

    if ($nl == 0) {
      $worksheet->write(0,0,'Sample',$centerAlignHeader2);
      for (my $i=1;$i<@array;$i++) {
        $worksheet->write(0,$i,$array[$i],$centerAlignHeader2);
      }
    }
    else {
      $worksheet->write($nl,0,$array[0]);
      for (my $i=1;$i<@array;$i++) {
        $worksheet->write_number($nl,$i,$array[$i]);
      }
    }
    $nl ++;
  }
  close(IN);
  $worksheet->freeze_panes(1,1);

  $workbook->close();

  exit;
}

my @betaMatrix = ('weighted_unifrac','unweighted_unifrac');
if ($pwd =~ /ITS/) {
  @betaMatrix = ('euclidean','pearson');
}

if (fork() == 0) {
  my @sampleOrder;
  my %collect;
  foreach my $matrix (@betaMatrix) {
    open(IN,'<'.$output.$matrix.'_'.$outputFileName);
    chomp($_ = <IN>);
    @sampleOrder = split("\t",$_);
    my @header = split("\t",$_,-1);
    while (<IN>) {
      chomp;
      if ($_ eq '') {
        next;
      }

      my @array = split("\t",$_);

      for (my $i=1;$i<@array;$i++) {
        $collect{$matrix}{$array[0]}{$header[$i]} = $array[$i];
      }
    }
    close(IN);
  }

  my $workbook = Spreadsheet::WriteExcel->new($output.'betaDiversity.xls');
  my $centerAlignHeader1 = $workbook->add_format(valign=>'vcenter',align=>'center',bold=>1);
  my $centerAlignHeader1Rotate = $workbook->add_format(valign=>'vcenter',align=>'center',bold=>1,rotation=>90);
  my $centerAlignHeader2 = $workbook->add_format(align=>'center',bold=>1);
  my $centerAlign = $workbook->add_format(valign=>'vcenter',align=>'center');

  my @formatWeighted;
  foreach my $colorIndex (0..10) {
    $workbook->set_custom_color(($colorIndex+10),255,25.5*(10-$colorIndex),25.5*(10-$colorIndex)); # 10 .. 20 for weighted
    $formatWeighted[$colorIndex] = $workbook->add_format(bg_color => ($colorIndex+10));
  }
  $formatWeighted[0] = $workbook->add_format(bg_color => 10);
  $formatWeighted[10] = $workbook->add_format(bg_color => 20);

  my @formatUnweighted;
  foreach my $colorIndex (0..10) {
    $workbook->set_custom_color(($colorIndex+30),25.5*(10-$colorIndex),255,25.5*(10-$colorIndex)); # 30 .. 40 for unweighted
    $formatUnweighted[$colorIndex] = $workbook->add_format(bg_color => ($colorIndex+30));
  }
  $formatUnweighted[0] = $workbook->add_format(bg_color => 30);
  $formatUnweighted[10] = $workbook->add_format(bg_color => 40);

  my $worksheet = $workbook->add_worksheet('beta');

  for (my $row=0;$row<@sampleOrder;$row++) {
    for (my $col=0;$col<@sampleOrder;$col++) {
      if (($row == 0) && ($col == 0)) {
        next;
      }
      elsif ($row == 0) {
        $worksheet->write($row,$col,$sampleOrder[$col],$centerAlignHeader2);
      }
      elsif ($col == 0) {
        $worksheet->write($row,$col,$sampleOrder[$row],$centerAlignHeader2);
      }
      else {
        if ($row > $col) { # bottom-left
          my $val = $collect{$betaMatrix[0]}{$sampleOrder[$row]}{$sampleOrder[$col]};
          my $colorIndex = 10*$val;

          $worksheet->write_number($row,$col,$val,$formatWeighted[$colorIndex]);
        }
        elsif ($row < $col) { #top-right
          my $val = $collect{$betaMatrix[1]}{$sampleOrder[$row]}{$sampleOrder[$col]};
          my $colorIndex = 10*$val;

          $worksheet->write_number($row,$col,$val,$formatWeighted[$colorIndex]);
        }
        elsif ($row == $col) {
          my $val = $collect{$betaMatrix[0]}{$sampleOrder[$row]}{$sampleOrder[$col]};
          my $colorIndex = 10*$val;

          #$worksheet->write($row,$col,$val,$formatWeighted[$colorIndex]);
          $worksheet->write($row,$col,'--');
        }
      }
    }
  }

  my $nl = @sampleOrder;
  $worksheet->merge_range($nl,1,$nl,($nl-1),$betaMatrix[0],$centerAlignHeader1);
  $worksheet->merge_range(1,$nl,($nl-1),$nl,$betaMatrix[1],$centerAlignHeader1Rotate);


  #$centerAlignHeader2

  $worksheet->freeze_panes(1,1);
  $workbook->close();

  exit;
}

wait for 0 .. 2;

