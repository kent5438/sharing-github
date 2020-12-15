#!/usr/bin/perl
use strict;
use warnings;
use File::Path qw(make_path);

my $input_taxa = defined $ARGV[0]? $ARGV[0] : '';
my $output_folder = defined $ARGV[1]? $ARGV[1] : '';

die "perl $0 [input taxa txt] [output folder]\n" if(!-f $input_taxa || $output_folder eq "");

$output_folder =~ s/\///g;
make_path($output_folder);

# read taxa file
my $defaultHeader = 1;
my @taxa = ();
my %sampleName = ();
my %sampleData = ();
my $leastColumn = 0;
open taxaFile, $input_taxa;
	while(<taxaFile>){
		#if($defaultHeader){$defaultHeader = 0; next;} # skip first line
		chomp;
		my @column = split(/\t/, $_);
		if($_ =~ /^Taxon/){ # header line
			for my $i (1..$#column){
				$sampleName{$i} = $column[$i];
				
				# init %$sampleData
				my @data = ();
				$sampleData{$i} = \@data
			}
			$leastColumn = $#column;
		}else{
			next if($#column < $leastColumn);
            
			$column[0] =~ s/;/\t/g; # translate ; to tab
			push @taxa, $column[0];
		    for my $i (1..$#column){
				push @{$sampleData{$i}}, $column[$i];
			}
		}
	}
close taxaFile;



# Output - one sample one output
foreach my $index (sort{$a<=>$b} keys %sampleName){
    open(outputFile, ">$output_folder/$sampleName{$index}.txt");
    my @data = @{$sampleData{$index}};
    for my $i(0..$#data){
        print outputFile $data[$i] ."\t". $taxa[$i]."\n";
    }
    close outputFile;
}
