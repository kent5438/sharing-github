#!/usr/bin/perl -w
#
# Process VEP output (VCF)
# /export/EC1680U/Daniel_WES/GNND-Zheng160127/DNA15-064/output/vqsr
# variant_effect_predictor.pl -i DNA15-064.filtered.vcf --cache --everything --dir /export/arrayPRO1/commonSoftwares/ensembl-tools-release-81/.vep/ --cache_version 81 --pick --vcf --fork 30 --ASSEMBLY GRCh38
# ./vcf2txt.v2.pl variant_effect_output.txt > DNA15-064.variant_effect_output.txt

my $fn = shift; # ouput file of VEP
my $header = ""; # VCF header
my $INFO; # fileds of VEP INFO
my @samples;
my @SamplesField;
my @fields;
my $fields; # fields name
my $sampledata; # hash structure with allele depth, allele frequency of each sample
my $finalfield = "";
my $sentinel = 1; # 哨兵
my $record; # one record of VEP output
my $string;
my $sample = "";
open (IN, "<$fn") or die "Can't open $fn file";
while (<IN>) {
  $ln = $_;
  chomp($ln);
  if ($ln=~/#/) { # VCF header
    if ($ln=~/#CHROM/) {
      # 取出 merged vcf 樣品名稱
      my @tab = split(/\t/, $ln);
      for (my $i=9; $i < $#tab + 1; $i++) {
        push (@samples, $tab[$i]);
        $string = $tab[$i]." frequency";
        push (@SamplesField, $string);
        $string = $tab[$i]." depth";
        push (@SamplesField, $string);
      }
      $sampledata = build(\@samples); # arrary reference
      next;
    }
    $header = $header.$ln."\n";
    if ($ln=~/##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: (.*)\">/) {
      $INFO = $1;
      @fields = split(/\|/, $INFO);
      $fields = processINFO(\@fields);
    }
  } else { # VCF record
    if ($sentinel ==1) {
      #$finalfield = "#CHROM\tPOS\tID\tREF\tALT\t".join("\t", @SamplesField)."\tSET\tQUAL\tFILTER\t".$fields;
      $finalfield = "#CHROM\tPOS\tID\tREF\tALT\t".join("\t", @SamplesField)."\tQUAL\tFILTER\t".$fields;
      print "$header$finalfield";
      $sentinel = $sentinel + 1;
    }
    $record = processRecord($ln, \@fields, \@samples, $sampledata);
    print "$record\n";
  }
}
close IN;

sub processINFO {
  # This process maybe useless when process New VEP output, because New VEP output will not have Allele base info
  # But this may be useful for Plugin (Condel) or other filed output
  my $fields_ref = shift;
  my @fields = @$fields_ref;
  my $fields = "";
  my $idx;
  for ($idx=1; $idx < scalar @fields; $idx++) {
    if ($fields[$idx] eq "Condel") {
      $fields = $fields."Condel_Desc\tCondel_value\t";
    } elsif ($fields[$idx] eq "SIFT") {
      $fields = $fields."SIFT_Desc\tSIFT_value\t";
    } elsif ($fields[$idx] eq "PolyPhen") {
      $fields = $fields."PolyPhen_Desc\tPolyPhen_value\t";
    } elsif ($fields[$idx] eq "GMAF") {
      $fields = $fields."GMAF_Allele\tGMAF_Freq\t";
    } elsif ($fields[$idx] eq "AFR_MAF") {
      $fields = $fields."AFR_MAF_Allele\tAFR_MAF_Freq\t";
    } elsif ($fields[$idx] eq "AMR_MAF") {
      $fields = $fields."AMR_MAF_Allele\tAMR_MAF_Freq\t";
    } elsif ($fields[$idx] eq "ASN_MAF") {
      $fields = $fields."ASN_MAF_Allele\tASN_MAF_Freq\t";
    } elsif ($fields[$idx] eq "EAS_MAF") {
      $fields = $fields."EAS_MAF_Allele\tEAS_MAF_Freq\t";
    } elsif ($fields[$idx] eq "EUR_MAF") {
      $fields = $fields."EUR_MAF_Allele\tEUR_MAF_Freq\t";
    } elsif ($fields[$idx] eq "SAS_MAF") {
      $fields = $fields."SAS_MAF_Allele\tSAS_MAF_Freq\t";
    } elsif ($fields[$idx] eq "AA_MAF") {
      $fields = $fields."AA_MAF_Allele\tAA_MAF_Freq\t";
    } elsif ($fields[$idx] eq "EA_MAF") {
      $fields = $fields."EA_MAF_Allele\tEA_MAF_Freq\t";
    } else {
      $fields = $fields.$fields[$idx]."\t";
    }
  }
  chop($fields);
  $fields = $fields."\n";
  return $fields;
}

sub build {
  my $samples = shift; # arrary reference
  my $data;
  my $vcf;
  foreach my $sample (@$samples) {
    $vcf = $sample.".vcf";
    open (VCF, '<', $vcf);
    while (<VCF>) {
      my $ln = $_;
      chomp($ln);
      if ($ln=~/#/) {next;}

      my @tab = split(/\t/, $ln);
      $chr=$tab[0];
      $pos=$tab[1];
      $reb=$tab[3]; # reference base
      $alb=$tab[4]; # allele base
      $format=$tab[-1];
      my @colon = split(/:/, $format);                # FORMAT GT:AD:DP:GQ:PL
      (my $refcount, my $altcount) = split(/,/, $colon[1]);
      my $depth = $colon[2];
      #print "$refcount\t$altcount\n";
      if ($refcount+$altcount == 0) {
        # 有些 variats Format 的 AD 會是 0,0, 原因未知
        #
        #chr13   93096604        .       T       TTCTATGAACCATGGGTTCAGGCGCATGCTCCTCCTCTCTGTGGACTATGGGCTCAGGCGCATGCTCCTCCTCTCTGCGGACTATGGGCTCAGGTGCATGCTCTTC      1544.73 .       AC=2;AF=1.00;AN=2;DP=106;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=55.79;SOR=0.693       GT:AD:DP:GQ:PL  1/1:0,0:0:65:1582,65,0
        #print "$chr\t$pos\n";
        #exit;
        $data->{$sample}->{$chr}->{$pos}->{reb} = $reb;
        $data->{$sample}->{$chr}->{$pos}->{alb} = $alb;
        $data->{$sample}->{$chr}->{$pos}->{dep} = 'NA';
        $data->{$sample}->{$chr}->{$pos}->{freq} = 'NA';
        next;
      }
      my $frequency = sprintf("%.3f",$altcount/($refcount+$altcount));
      $data->{$sample}->{$chr}->{$pos}->{reb} = $reb;
      $data->{$sample}->{$chr}->{$pos}->{alb} = $alb;
      $data->{$sample}->{$chr}->{$pos}->{dep} = $depth;
      $data->{$sample}->{$chr}->{$pos}->{freq} = $frequency;
      #print "$sample\t$chr\t$pos\t$depth\t$frequency\n";
    }
    close VCF;
  }
  return $data;
}

sub processRecord {
  my $ln = shift; # one record of VEP output
  my $fields_ref = shift;
  my $samples_ref = shift; # arrary reference
  my $sampledata = shift;
  my @fields = @$fields_ref;
  my @samples = @$samples_ref;
  my @tab  = split(/\t/, $ln);
  my $chr  = $tab[0];                 # chromosome
  my $pos  = $tab[1];                 # position
  my $id   = $tab[2];                 # ID
  my $ref  = $tab[3];                 # Reference Allele
  my $alt  = $tab[4];                 # Alternative Allele
  my $qual = $tab[5];                 # Quality
  my $filt = $tab[6];                 # filter
  my $info = $tab[7];                 # Info
  my $fielddata;
  #$info=~/.*(set=.*);.*/;
  #my $set = $1;
  my @semicolon = split(/;/, $info);

  my $tmp   = $semicolon[-1];
  $tmp =~/CSQ=(.*)/;
  my $csq = $1;
  $csq =~s/\|/ \|/g;     # Very Important 用 space (空格) 塞入沒有資料的欄位
  my @pipe = split(/\|/, $csq);
  for (my $idx=1; $idx < scalar @fields; $idx++) { # escape Allele information
    my $tmp_value = $pipe[$idx];
    unless (defined $pipe[$idx]) {
      $tmp_value = " ";
    }
    $fielddata->{$fields[$idx]} = $tmp_value; # Each Field has its value
  }

  my @SamplesInfo = ();
  foreach my $sample (@samples) {
    if (defined $sampledata->{$sample}->{$chr}->{$pos}) {
      $reb       = $sampledata->{$sample}->{$chr}->{$pos}->{reb};
      $alb       = $sampledata->{$sample}->{$chr}->{$pos}->{alb};
      $depth     = $sampledata->{$sample}->{$chr}->{$pos}->{dep};
      $frequency = $sampledata->{$sample}->{$chr}->{$pos}->{freq};
      #print "$sample\t$reb\t$alb\t$depth\t$frequency\t$ref\t$alt\n";
      if (($ref eq $reb) && ($alt eq $alb)) {
        push (@SamplesInfo, $frequency); # frequency
        push (@SamplesInfo, $depth); # depth
      } else {
        # 相同位置, 不同 allele base
        push (@SamplesInfo, 'NA'); # frequency
        push (@SamplesInfo, 'NA'); # depth
      }
    } else {
      push (@SamplesInfo, 'NA'); # frequency
      push (@SamplesInfo, 'NA'); # depth
    }
  }
  #"#CHROM\tPOS\tID\tREF\tALT\t".join("\t", @SamplesField)."\tSET\tQUAL\tFILTER\t".$fields;
  #my $output = "$chr\t$pos\t$id\t$ref\t$alt\t".join("\t", @SamplesInfo)."\t$set\t$qual\t$filt\t";
  my $output = "$chr\t$pos\t$id\t$ref\t$alt\t".join("\t", @SamplesInfo)."\t$qual\t$filt\t";
  for ($idx=1; $idx < scalar @fields; $idx++) {
    my $tmp_output = $fielddata->{$fields[$idx]};
    if ($fields[$idx] eq "Condel") {
      $newoutput = processvalue($tmp_output);
    } elsif ($fields[$idx] eq "SIFT") {
      $newoutput = processvalue($tmp_output);
    } elsif ($fields[$idx] eq "PolyPhen") {
      $newoutput = processvalue($tmp_output);
    } elsif ($fields[$idx] eq "GMAF") {
      $newoutput = processvalue($tmp_output);
    } elsif ($fields[$idx] eq "AFR_MAF") {
      $newoutput = processvalue($tmp_output);
    } elsif ($fields[$idx] eq "AMR_MAF") {
      $newoutput = processvalue($tmp_output);
    } elsif ($fields[$idx] eq "ASN_MAF") {
      $newoutput = processvalue($tmp_output);
    } elsif ($fields[$idx] eq "EAS_MAF") {
      $newoutput = processvalue($tmp_output);
    } elsif ($fields[$idx] eq "EUR_MAF") {
      $newoutput = processvalue($tmp_output);
    } elsif ($fields[$idx] eq "SAS_MAF") {
      $newoutput = processvalue($tmp_output);
    } elsif ($fields[$idx] eq "AA_MAF") {
      $newoutput = processvalue($tmp_output);
    } elsif ($fields[$idx] eq "EA_MAF") {
      $newoutput = processvalue($tmp_output);
    } else {
      $newoutput = $tmp_output;
    }
    $output = $output.$newoutput."\t";
  }
  chop($output);
  return $output;
}
sub processvalue {
  my $value = shift;
  my $newvalue;
  if ($value=~/(\D+)\((.*)\)/) {
    $newvalue = "$1\t$2\t";
  } elsif ($value=~/(\D)\:(.*)/) {
    $newvalue = "$1\t$2\t";
  } else {
    $newvalue = "\t\t";
  }
  return $newvalue;
}
1;
