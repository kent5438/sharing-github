#!/usr/bin/perl -w
#
# Copyright belong to Chen Yuan Liu 
# SCRIPT STATEMENT
# SCRIPT STATEMENT
# Author: Chen-Yuan Liu
# Date: 03-04, 2016
# Data structure
#  
# Function:
# ./heatmap.pl --comps=output/ebseq/AM1-D7,AM3-D7-AM1-2D-1.genMat.result.1#output/ebseq/AM1-D7,AM3-D7-AM3-2D-1.genMat.result.1
#

$main = new heatmap;
$main-> parse_command_line;
$main-> geneheatmap;
$main-> export;

package heatmap;
use strict;
use warnings;
use Getopt::Long;
use threads;
use vars qw($app %parameters %data @program);

sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;
  $self->init;
  return $self;
}
sub init {
  my $self = shift;
}
sub parse_command_line {
  my $self = shift;
  our $parameters;
  my $comps = "";  # 比對方案
  my $output = ""; # output folder

  if(scalar(@ARGV) < 2){
    print usage();
	exit(1);
  } 
  GetOptions(
    'comps=s'    => \$comps,	
	'output=s'    => \$output,
  );
  $parameters->{comps} = $comps;
  $parameters->{output} = $output;
}
sub usage {
  my $self = shift;
  my $msg = <<_EOUSAGE_;
# REQUIRED:
# comps  <string>  :comparsion set
# output <string>  :output folder
_EOUSAGE_
;
  return $msg;
}
sub geneheatmap {
  my $self = shift;
  our $parameters;
  my $comps = $parameters->{comps};
  my @comma;
  my $comma;
  my $fn;
  @comma = split(/\#/, $comps);
  foreach $comma (@comma) {
	$self -> build($comma);
  }
}
sub build {   
  our $data;  
  our @comparsions;
  
  my $self = shift;
  my $fn   = shift;
  my $comparsion;
  my $ln;
  my @tab;
  my @slash;
  my $gene;
  my $ppde;
  my $postfc;
  my $program;
  @slash = split(/\//, $fn);
  $comparsion = $slash[-1];  
  push (@comparsions, $comparsion);
  open (IN, "<$fn") or die "Can't open $fn file";
  while (<IN>) {
    $ln = $_;
	chomp($ln);
	if ($ln=~/PPEE/) {
	  next;
	}
	@tab = split(/\t/, $ln);
	$gene = $tab[0];
	$ppde = $tab[2];
	$postfc = $tab[3];	
#	if ($ppde > 0.95) {
#	  if (($postfc < 0.5) or ($postfc > 2 )) {	    
	    $data->{$gene}-> {$comparsion}->{FC} = $postfc;
#	  }
#	}
  }
  close IN;  
}
sub export {  
  our @comparsions;
  our $data;
  our $parameters;
  my $self = shift;
  my $comparsion;
  my $gene;
  my $fc;
  my $outdir = $parameters->{output};
  my $output = "$outdir/comparsion_union.txt";
  my $check;
  open (OUT, ">$output") or die "Can't create $output file\n";
  print OUT "Gene\t";
  foreach $comparsion (@comparsions) {
    print OUT "$comparsion\t";
  }
  print OUT "\n";
  
  foreach $gene (sort keys %{$data}) {
    $check = 1;
	
	foreach $comparsion (@comparsions) {
	  unless (defined $data->{$gene}->{$comparsion}) { 
	    $check = 0;
	  }
	}	
	if ($check) {
	  print OUT "$gene\t";
	  foreach $comparsion (@comparsions) {
	    $fc = $data->{$gene}-> {$comparsion}->{FC};
		print OUT "$fc\t";
	  }
	  print OUT "\n";
	}	
  }  
  close OUT;
}
1;
