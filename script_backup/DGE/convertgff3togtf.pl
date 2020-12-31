#!/usr/bin/perl -w

$gff3 = shift;
$gtf = $gff3;
$gtf=~s/gff3/gtf/;

$gene = "";
$transcript_id = "";
open (OUT, '>', $gtf) or die "Can't create $gtf file\n"; 
open (IN, '<',  $gff3) or die "Can't read $gff3 file\n";
while (<IN>) {
  $ln = $_;
  chomp($ln);
  @tab = split(/\t/, $ln);
  if ($ln=~/^#/) {next;}
  if ($tab[2] eq "gene") {
    $tab[8]=~/ID=(.*)/;
    @tmp = split(";", $1);
    $gene = $tmp[0];
  } 
  #if ($tab[2] eq "transcript") {
  if ($tab[2] eq "mRNA") {
    $tab[8]=~/ID=(.*);Parent=$gene/;  
    $transcript_id = $1;
  } 
  if (($tab[2] eq "CDS") or ($tab[2] eq "exon")) {
    print OUT  "$tab[0]\t$tab[1]\t$tab[2]\t$tab[3]\t$tab[4]\t0.000000\t$tab[6]\t.\tgene_id \"$gene\"; transcript_id \"$transcript_id\";\n";
  }
}
close IN;
close OUT;
