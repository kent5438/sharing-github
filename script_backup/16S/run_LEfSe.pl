#!/usr/bin/perl
use strict;
use warnings;
use File::Path qw(make_path);

my $input_taxa	= defined $ARGV[0]? $ARGV[0] : '';
my $mapping_file= defined $ARGV[1]? $ARGV[1] : '';
my $group_title	= defined $ARGV[2]? $ARGV[2] : 'Description';
my $output_folder=defined $ARGV[3]? $ARGV[3] : 'LEfSe';

die "perl $0 [input taxa txt] [mapping txt] ([group title: Description]) ([output folder: LEfSe])\n" if(!-f $input_taxa || !-f $mapping_file);
make_path($output_folder);


# read mapping file for group
my %group = ();
my %group_count = ();
my $idIndex = 0;
my $groupIndex = 0;
open mapFile, $mapping_file;
	while(<mapFile>){
		chomp;
		my @column = split(/\t/, $_);
		next if($#column < 3);
		if($_ =~ /^#/){
			for my $i(0..$#column){
				$idIndex = $i if($column[$i] eq '#SampleID');
				$groupIndex = $i if($column[$i] eq $group_title);
			}
		}else{
			$column[$groupIndex] =~ s/\s/_/g; # replace space
			$group{$column[$idIndex]} = $column[$groupIndex];
            $group_count{$column[$groupIndex]}++;
		}
	}
close mapFile;

die "LEfSe need two groups at least.\n" if(keys %group_count < 2);

# read taxa file
my $least = (keys %group);
my %taxa;
my @sample = ();
open taxaFile, $input_taxa;
	while(<taxaFile>){
		chomp;
		my @column = split(/\t/, $_); # OTU, samples..
		next if($#column < $least);
		
		#if($_ =~ /^#/){
		if($_ =~ /^Taxon/){
			@sample = @column[1..$#column];
		}else{
			$column[0] =~ s/\s//g; # replace space
			$column[0] =~ s/;/\|/g; # replace spliter
			next if($column[0] =~ /[kpcofgs]__$/); # skip unknow case
			
			$column[0] =~ s/^[kpcofgs]__//g;    # replace k__, p__, c__, o__, f__, g__, s__
			$column[0] =~ s/\|[kpcofgs]__/\|/g; # replace k__, p__, c__, o__, f__, g__, s__
			
			# merge if exist
			if(defined($taxa{$column[0]})){
				my @cache = split("\t", $taxa{$column[0]});
				for my $i(1..$#column){
					$column[$i] += $cache[$i-1];
				}
			}
			$taxa{$column[0]} = join("\t", @column[1..$#column]);
		}
	}
close taxaFile;

# output
open(outFile, ">$output_folder/LEfSe_input.txt");
	print outFile join("\t", 'SampleID', @sample) ."\n"; # print sample 
	
	my @subTitle = ();
	for my $i(0..$#sample){
		push @subTitle, $group{$sample[$i]};
	}
	print outFile join("\t", $group_title, @subTitle) ."\n"; # print subset title
	
	for my $otu (keys %taxa){
		print outFile join("\t", $otu, $taxa{$otu}) ."\n"; # print taxa
	}
close outFile;



# run LEfSe
`format_input.py $output_folder/LEfSe_input.txt $output_folder/LEfSe_formated.txt -u 1 -c 2 -o 1000000`;
my $log = `run_lefse.py $output_folder/LEfSe_formated.txt $output_folder/LEfSe_output.txt`;
die "No features with significant differences between the two classes\n" if($log =~ /No features with significant differences between the two classes/);
`plot_res.py $output_folder/LEfSe_output.txt $output_folder/LEfSe.png --width 15 --feature_font_size 10`;
`plot_cladogram.py $output_folder/LEfSe_output.txt  $output_folder/LEfSe_cladogram.png --format png --dpi 150  --left_space_prop 0`;
