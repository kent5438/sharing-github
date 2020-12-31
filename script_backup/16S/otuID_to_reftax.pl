#!/usr/bin/env perl
use strict;
# author: jsyf, 2018-05-10

=h
my @taxLevel = qw / k p c o f g s / ;

my %silvaTaxToId;
open(IN,'<./raw_taxonomy.txt');
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  my ($id,$tax) = split("\t",$_);

  #my @allTax = split(';',$tax);

  #for (my $i=0;$i<@allTax;$i++) {
  #  if ($allTax[$i] =~ /^D_${i}__/) {
  #    $allTax[$i] =~ s/^D_${i}__/$taxLevel[$i]__/;
  #  }
  #  #else { # Ambiguous_taxa
  #  #  print $id."\t".$tax."\n";
  #  #}
  #}

  #print join("   ",$id,join(' , ',@allTax))."\n";
  if (!defined($silvaTaxToId{$tax})) {
    $silvaTaxToId{$tax} = $id;
  }
}
close(IN);
=cut

foreach my $file ('trimmed.unique.opti_mcc.0.03.rep.fasta','results.blast6','otu_table.otuID.filter.dense.txt') {
  if (!-e $file) {
    print $file.' not found!'."\n";
    exit;
  }
}

my %otuToSeqId;
open(IN,'<trimmed.unique.opti_mcc.0.03.rep.fasta');
while (my $name = <IN>){
  chomp($name);
  chomp(my $seq = <IN>);

  if ($name =~ /^>(\S+)\s+([^\|]+)\|/) {
    my ($seq,$otu) = ($1,$2);

    #print $seq."\t".$otu."\n";
    $otuToSeqId{$otu} = $seq;
  }
  else {
    next;
  }
}
close(IN);

my %seqToSilva;
open(IN,'<results.blast6');
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  my @array = split("\t",$_);

  $seqToSilva{$array[0]} = $array[1];
}
close(IN);

open(OUT,'>otu_table.taxID.filter.txt');
open(FAIL,'>otu_table.taxID.filter.txt_failed');
open(IN,'<otu_table.otuID.filter.txt');
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  if ($_ =~ /^#/) {
    print OUT $_."\n";
    next;
  }

  my @array = split("\t",$_);
  my $otu = shift(@array);

  if (!defined($otuToSeqId{$otu}) || !defined($seqToSilva{$otuToSeqId{$otu}})) {
    print FAIL join("\t",$otu,@array)."\n";
  }
  else {
    my $silva = $seqToSilva{$otuToSeqId{$otu}};

    print OUT join("\t",$silva,@array)."\n";
  }
}
close(IN);
close(FAIL);
close(OUT);

use JSON;
my $str = '';

open(IN,'<otu_table.otuID.filter.biom');
while (<IN>) {
  chomp;
  $str .= $_;
}
close(IN);

open(FAIL,'>otu_table.taxID.filter.biom_failed');
my $json_ref = decode_json($str);
for (my $i=0;$i<@{$json_ref->{'rows'}};$i++) {
  my $otu = $json_ref->{'rows'}[$i]{'id'};

  if (!defined($otuToSeqId{$otu}) || !defined($seqToSilva{$otuToSeqId{$otu}})) {
    print FAIL $otu."\n";
  }
  else {
    $json_ref->{'rows'}[$i]{'id'} = $seqToSilva{$otuToSeqId{$otu}};
  }
}
close(FAIL);

open(OUT,'>otu_table.taxID.filter.dense.biom');
my $output = encode_json($json_ref);
print OUT $output."\n";
close(OUT);

=h
==> ./majority_taxonomy_7_levels.txt <==
KF494428.1.1396	D_0__Bacteria;D_1__Epsilonbacteraeota;D_2__Campylobacteria;D_3__Campylobacterales;D_4__Thiovulaceae;D_5__Sulfuricurvum;D_6__Sulfuricurvum sp. EW1
AF506248.1.1375	D_0__Bacteria;D_1__Cyanobacteria;D_2__Oxyphotobacteria;D_3__Nostocales;D_4__Nostocaceae;D_5__Nostoc PCC-73102;D_6__Nostoc sp. 'Nephroma expallidum cyanobiont 23'
EF603722.1.1487	D_0__Bacteria;D_1__Bacteroidetes;D_2__Bacteroidia;D_3__Bacteroidales;D_4__Muribaculaceae;D_5__uncultured bacterium;D_6__uncultured bacterium

==> ./otu_table.otuID.filter.txt <==
# Constructed from biom file
#OTU ID	AY457674	AY457675	taxonomy
Otu01	0.0	0.0	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__Faecalibacterium; g__Faecalibacterium_unclassified

>AY457780	Otu01|20|AY457695-AY457732-AY457748-AY457754-AY457755-AY457758-AY457763-AY457767-AY457780-AY457785-AY457793-AY457800-AY457809-AY457818-AY457826-AY457847-AY457849-AY457860-AY457871-AY457881
CCCTTAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGCGCCTAACACATGCAAGTCGAACGAGCGAGAGAGAGCTTGCTTTCTCGAGCGAGTGGCGAACGGGTGAGTAACGCGTGAGGAACCTGCCTCAAAGAGGGGGACAACAGTTGGAAACGACTGCTAATACCGCATAAGCCCACGGGTCGGCATCGACCAGAGGGAAAAGGAGCAATCCGCTTTGAGATGGCCTCGCGTCCGATTAGCTAGTTGGTGAGGTAACGGCCCACCAAGGCAACGATCGGTAGCCGGACTGAGAGGTTGAACGGCCACATTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGCGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGGAGGAAGAAGGTCTTCGGATTGTAAACTCCTGTTGTGAGGAAGATAAAGACGGTACTCAACAAGG

==> results.blast6 <==
AY457740	CDYU01016198.40941.42445	32.3	495	3	3	1	497	1	1505	-1	0
AY457686	GQ493387.1.1374	24.5	340	2	2	1	338	1	1374	-1	0
AY457848	EF404814.1.1488	32.6	490	3	0	1	495	1	1487	-1	0

            {"id":"Otu03", "metadata":{"taxonomy":["k__Bacteria", "p__Bacteroidetes", "c__Bacteroidia", "o__Bacteroidales", "f__Bacteroidaceae", "g__Bacteroides", "s__Bacteroides eggerthii DSM 20697"]}},
            {"id":"Otu04", "metadata":{"taxonomy":["k__Bacteria", "p__Bacteroidetes", "c__Bacteroidia", "o__Bacteroidales", "f__Tannerellaceae", "g__Parabacteroides"]}},
            {"id":"Otu05", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Betaproteobacteriales", "f__Burkholderiaceae", "g__Sutterella"]}},
            {"id":"Otu06", "metadata":{"taxonomy":["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Ruminococcaceae", "g__Faecalibacterium", "g__Faecalibacterium_unclassified"]}},
            {"id":"Otu07", "metadata":{"taxonomy":["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Lachnospiraceae", "g__Lachnospira"]}},
=cut
