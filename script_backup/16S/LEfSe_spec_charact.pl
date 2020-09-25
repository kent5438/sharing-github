#! /usr/bin/env perl

use strict;

my $docker = "sudo docker run -it --rm  -v /etc/localtime:/etc/localtime:ro -w `pwd` -v /export/EC1680U/:/export/EC1680U/ -v /export/EC1680U/DataBase/16S:/references -v `pwd`:`pwd` kentchendocker/g16s_pipeline:latest";

if(! -e "../D.1_taxa_summary/otu_table.taxID.filter.dense.dedup_L6.txt.bak"){
	system("cp ../D.1_taxa_summary/otu_table.taxID.filter.dense.dedup_L6.txt ../D.1_taxa_summary/otu_table.taxID.filter.dense.dedup_L6.txt.bak");
}

open(OUT,'>../D.1_taxa_summary/otu_table.taxID.filter.dense.dedup_L6.txt');
open(IN,'<../D.1_taxa_summary/otu_table.taxID.filter.dense.dedup_L6.txt.bak');
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  my @array = split("\t",$_);

  if ($. == 1) {
    print OUT join("\t",@array)."\n";
    next;
  }

  my $tax = shift(@array);

  $tax =~ s/[^\w-\t\.(?!;$)]+/_/g;

  print OUT join("\t",$tax,@array)."\n";
}
close(IN);
close(OUT);

if(-e "../../data/Mapping_1.txt"){
	system("$docker /export/EC1680U/perl/bin/16S/run_LEfSe.pl ../D.1_taxa_summary/otu_table.taxID.filter.dense.dedup_L6.txt ../../data/Mapping_1.txt Description .");
} else {
	system("$docker /export/EC1680U/perl/bin/16S/run_LEfSe.pl ../D.1_taxa_summary/otu_table.taxID.filter.dense.dedup_L6.txt ../../data/Mapping.txt Description .");
}
