#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);

my $input_otu	= defined $ARGV[0]? $ARGV[0] : '';
my $mapping_file= defined $ARGV[1]? $ARGV[1] : '';
my $output_file	= defined $ARGV[2]? $ARGV[2] : 'OTU_venn.png';
my $group_title	= defined $ARGV[3]? $ARGV[3] : 'Description';

die "perl $0 [input otu txt] [mapping txt] ([output: OTU_venn.png]) ([group title: Description])\n" if(!-f $input_otu || !-f $mapping_file);
make_path(dirname($output_file));

# set color
our @color = ('#d85836', '#fac55e', '#79cada', '#bbd874', '#7e6bae');
our $alpha = 0.35;
our $catAdjust = 0;

# read mapping file for group
my %group = ();
my %otu = ();
my $idIndex = 0;
my $groupIndex = 0;
open mapFile, $mapping_file;
	while(<mapFile>){
		chomp;
		my @column = split(/\t/, $_);
		next if($#column < 3);
		if($_ =~ /^#(.+)/){
			@column = split(/\t/, $1);
			for my $i(0..$#column){
				$idIndex = $i if($column[$i] eq 'SampleID');
				$groupIndex = $i if($column[$i] eq $group_title);
			}
		}else{
			$group{$column[$idIndex]} = $column[$groupIndex];
			$catAdjust = (length($column[$groupIndex])>10 || $catAdjust)? 1 : 0;
			
			# init otu
			my %cache = ();
			$otu{$column[$groupIndex]} = \%cache;
		}
	}
close mapFile;

# read otu table
my $least = (keys %group);
my @sample = ();
open otuFile, $input_otu;
	while(<otuFile>){
		chomp;
		my @column = split(/\t/, $_); # OTU, samples.., taxa
		next if($#column < $least);
		
		if($_ =~ /^#/){
			@sample = @column[1..$#column-1];
		}else{
			for my $i(1..$#column-1){
				$otu{$group{$sample[$i-1]}}{$column[0]} = 1 if($column[$i]>0);
			}
		}
	}
close otuFile;

&DrawVenn(\%otu, $output_file);




#-- Function part --#

sub DrawVenn(){
	my %group = %{shift @_};
	my $fileName = shift @_;
		
	# setup function
	my $code = '';
	my $total = keys %group;
	if($total == 1){
		$code .= &OneVenn(\%group);
	}elsif($total == 2){
		$code .= &TwoVenn(\%group);
	}elsif($total == 3){
		$code .= &ThreeVenn(\%group);
	}elsif($total == 4){
		$code .= &FourVenn(\%group);
	}elsif($total == 5){
		$code .= &FiveVenn(\%group);
	}else{
		#die "There has $total group, total group should be 1 ~ 5\n";
                exit;
	}

	# create R script
	open(scriptFile, ">venn_cache.R");
		print scriptFile "library(VennDiagram)\n";
		print scriptFile "png(file='$fileName', width=1080, height=1080)\n";;
		print scriptFile "$code\n";
		print scriptFile "dev.off()";
	close scriptFile;
	
	`Rscript venn_cache.R`;
}

sub OneVenn(){
	my %data = %{shift @_};
	my @name = sort{$a cmp $b} keys %data;
	my $area = keys %{$data{$name[0]}};
	
	return "draw.single.venn(area=$area, category='$name[0]', fill = '$color[0]', alpha=$alpha)";
}

sub TwoVenn(){
	my %data = %{shift @_};	
	my @name = sort{$a cmp $b} keys %data;
	my @area1 = keys %{$data{$name[0]}};
	my @area2 = keys %{$data{$name[1]}};
	
	my $code = 'draw.pairwise.venn(';
	$code .=   'area1='. ($#area1+1);
	$code .= ', area2='. ($#area2+1);
	$code .= ', cross.area='. &CrossArea(\@area1, \@area2);
	$code .= ", category=c('$name[0]', '$name[1]')";
	$code .= ", fill=c('$color[0]', '$color[1]')";
	$code .= ", alpha=$alpha";
	$code .= ", cat.just=list(c(0.2,-7) , c(1,-7))" if($catAdjust);
	$code .= ')';
	
	return $code;
}

sub ThreeVenn(){
	my %data = %{shift @_};	
	my @name = sort{$a cmp $b} keys %data;
	my @area1 = keys %{$data{$name[0]}};
	my @area2 = keys %{$data{$name[1]}};
	my @area3 = keys %{$data{$name[2]}};
	
	my $code = 'draw.triple.venn(';
	$code .=   'area1='. ($#area1+1);
	$code .= ', area2='. ($#area2+1);
	$code .= ', area3='. ($#area3+1);
	$code .= ', n12='. &CrossArea(\@area1, \@area2);
	$code .= ', n13='. &CrossArea(\@area1, \@area3);	
	$code .= ', n23='. &CrossArea(\@area2, \@area3);
	$code .= ', n123='. &CrossArea(\@area1, \@area2, \@area3);
	$code .= ", category=c('$name[0]', '$name[1]', '$name[2]')";
	$code .= ", fill=c('$color[0]', '$color[1]', '$color[2]')";
	$code .= ", alpha=$alpha";
	$code .= ", cat.just=list(c(-0.1,6), c(1.1,7), c(0.5,0))" if($catAdjust);
	$code .= ')';
	
	return $code;
}

sub FourVenn(){
	my %data = %{shift @_};	
	my @name = sort{$a cmp $b} keys %data;
	my @area1 = keys %{$data{$name[0]}};
	my @area2 = keys %{$data{$name[1]}};
	my @area3 = keys %{$data{$name[2]}};
	my @area4 = keys %{$data{$name[3]}};

	my $code = 'draw.quad.venn(';
	$code .=   'area1='. ($#area1+1);
	$code .= ', area2='. ($#area2+1);
	$code .= ', area3='. ($#area3+1);
	$code .= ', area4='. ($#area4+1);
	$code .= ', n12='. &CrossArea(\@area1, \@area2);
	$code .= ', n13='. &CrossArea(\@area1, \@area3);
	$code .= ', n14='. &CrossArea(\@area1, \@area4);	
	$code .= ', n23='. &CrossArea(\@area2, \@area3);
	$code .= ', n24='. &CrossArea(\@area2, \@area4);
	$code .= ', n34='. &CrossArea(\@area3, \@area4);
	$code .= ', n123='. &CrossArea(\@area1, \@area2, \@area3);
	$code .= ', n124='. &CrossArea(\@area1, \@area2, \@area4);
	$code .= ', n134='. &CrossArea(\@area1, \@area3, \@area4);
	$code .= ', n234='. &CrossArea(\@area2, \@area3, \@area4);
	$code .= ', n1234='. &CrossArea(\@area1, \@area2, \@area3, \@area4);
	$code .= ", category=c('$name[0]', '$name[1]', '$name[2]', '$name[3]')";
	$code .= ", fill=c('$color[0]', '$color[1]', '$color[2]', '$color[3]')";
	$code .= ", alpha=$alpha";
	$code .= ", cat.just=list(c(0.25,2), c(0.8,2), c(0.5,1), c(0.5,3))" if($catAdjust);
	$code .= ')';
	
	return $code;
}

sub FiveVenn(){
	my %data = %{shift @_};	
	my @name = sort{$a cmp $b} keys %data;
	my @area1 = keys %{$data{$name[0]}};
	my @area2 = keys %{$data{$name[1]}};
	my @area3 = keys %{$data{$name[2]}};
	my @area4 = keys %{$data{$name[3]}};
	my @area5 = keys %{$data{$name[4]}};
		
	my $code = 'draw.quintuple.venn(';
	$code .=   'area1='. ($#area1+1);
	$code .= ', area2='. ($#area2+1);
	$code .= ', area3='. ($#area3+1);
	$code .= ', area4='. ($#area4+1);
	$code .= ', area5='. ($#area5+1);
	$code .= ', n12='. &CrossArea(\@area1, \@area2);
	$code .= ', n13='. &CrossArea(\@area1, \@area3);
	$code .= ', n14='. &CrossArea(\@area1, \@area4);
	$code .= ', n15='. &CrossArea(\@area1, \@area5);
	$code .= ', n23='. &CrossArea(\@area2, \@area3);
	$code .= ', n24='. &CrossArea(\@area2, \@area4);
	$code .= ', n25='. &CrossArea(\@area2, \@area5);
	$code .= ', n34='. &CrossArea(\@area3, \@area4);
	$code .= ', n35='. &CrossArea(\@area3, \@area5);
	$code .= ', n45='. &CrossArea(\@area4, \@area5);
	$code .= ', n123='. &CrossArea(\@area1, \@area2, \@area3);
	$code .= ', n124='. &CrossArea(\@area1, \@area2, \@area4);
	$code .= ', n125='. &CrossArea(\@area1, \@area2, \@area5);
	$code .= ', n134='. &CrossArea(\@area1, \@area3, \@area4);
	$code .= ', n135='. &CrossArea(\@area1, \@area3, \@area5);
	$code .= ', n145='. &CrossArea(\@area1, \@area4, \@area5);
	$code .= ', n234='. &CrossArea(\@area2, \@area3, \@area4);
	$code .= ', n235='. &CrossArea(\@area2, \@area3, \@area5);
	$code .= ', n245='. &CrossArea(\@area2, \@area4, \@area5);
	$code .= ', n345='. &CrossArea(\@area3, \@area4, \@area5);
	$code .= ', n1234='. &CrossArea(\@area1, \@area2, \@area3, \@area4);
	$code .= ', n1235='. &CrossArea(\@area1, \@area2, \@area3, \@area5);
	$code .= ', n1245='. &CrossArea(\@area1, \@area2, \@area4, \@area5);
	$code .= ', n1345='. &CrossArea(\@area1, \@area3, \@area4, \@area5);
	$code .= ', n2345='. &CrossArea(\@area2, \@area3, \@area4, \@area5);
	$code .= ', n12345='. &CrossArea(\@area1, \@area2, \@area3, \@area4, \@area5);
	$code .= ", category=c('$name[0]', '$name[1]', '$name[2]', '$name[3]', '$name[4]')";
	$code .= ", fill=c('$color[0]', '$color[1]', '$color[2]', '$color[3]', '$color[4]')";
	$code .= ", alpha=$alpha";
	$code .= ", cat.just=list(c(0.5,1.5) , c(0,-6) , c(0.3,0) , c(0.7,1) , c(1.1,-6))" if($catAdjust);
	$code .= ')';

	return $code;
}

sub CrossArea(){
	my $total = $#_;
	
	my %total = ();
	foreach(@_){
		foreach my $item (@{$_}){
			$total{$item}++;
		}
	}
	
	my $cross = 0;
	foreach my $name (keys %total){
		$cross++ if($total{$name}>$total);
	}

	return $cross;
}
