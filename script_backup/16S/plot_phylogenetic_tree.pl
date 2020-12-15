#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);

my  $otu_table    = defined $ARGV[0]? $ARGV[0] : '';
my  $align_fasta  = defined $ARGV[1]? $ARGV[1] : '';
our $genera_top   = defined $ARGV[2]? $ARGV[2] : 10;
our $output_tree  = defined $ARGV[3]? $ARGV[3] : 'PhylogeneticTree.png';

die "perl $0 [otu table] [aligned fasta] ([top genera: 10]) ([output tree:PhylogeneticTree.png])\n" if(!-f $otu_table || !-f $align_fasta || $genera_top !~ /^[0-9]+$/);
make_path(dirname($output_tree));

# set argv
our @colorSet = ('#C45436', '#E9BA5E', '#73B9C7', '#AFC971', '#72639C', '#4AB9A1', '#4B5FAA', '#CE906D', '#AE6FAC', '#323031');
# fix all color
for my $i($#colorSet..$genera_top-2){
	push @colorSet, 'black';
}


# statistics otus
my $totalOTU  = 0;
my %sumGenera = ();
my %countOTU  = ();
my %groupOTU  = ();
open taxaFile, $otu_table;
	while(<taxaFile>){
		chomp;
		my @column = split(/\t/, $_); # OTU id, samples.., OTU taxa
		next if($#column < 2);
		
		if($_ !~ /^#/){
			$column[$#column] =~ s/\s//g;
			next if($column[$#column] =~ /g__;s__$/);
			
			my $sum = eval join('+', @column[1..$#column-1]);
			my @name = split(/;/, $column[$#column]); # k,p,c,o,f,g,s
			$name[5] =~ s/^g__//g;
			
			$groupOTU{$column[0]} = $name[5];
			$countOTU{$column[0]} = $sum;
			$sumGenera{$name[5]} += $sum;
			$totalOTU += $sum;
		}
	}
close taxaFile;

# pick top genera
my $count = 0;
my %topGenera = ();
foreach my $genera (reverse sort {$sumGenera{$a} <=> $sumGenera{$b}} keys %sumGenera){
	last if($count >= $genera_top);
	$topGenera{$genera} = $count;
	$count++;
}

# create annotate file
open(annFile, ">annotate_tree.txt");
	# global annotate
	print annFile join("\t", 'clade_marker_size', 0)."\n";
	print annFile join("\t", 'title', 'OTUs phylogenetic relationship')."\n";
	print annFile join("\t", 'title_font_size', 30)."\n";
	print annFile join("\t", 'class_legend_font_size', 20)."\n";
	print annFile join("\t", 'class_legend_marker_size', 4)."\n";
	print annFile join("\t", 'annotation_background_alpha', 0.8)."\n";
	
	# group annotate
	foreach my $genera (sort {$topGenera{$a} <=> $topGenera{$b}} keys %topGenera){
    my $trim = ($genera =~ /^(.+)_unclassified$/)? "*$1" : $genera;
		print annFile join("\t", $trim, 'clade_marker_color', $colorSet[$topGenera{$genera}])."\n";
	}
	
	# annotate edge color
	foreach my $id (keys %groupOTU){
		if(defined $topGenera{$groupOTU{$id}}){
			my $color = $colorSet[$topGenera{$groupOTU{$id}}];
			print annFile join("\t", $id, 'annotation_background_color', $color)."\n";
			print annFile join("\t", $id, 'ring_color', 1, $color)."\n";
			print annFile join("\t", $id, 'ring_height', 1, $countOTU{$id}/$totalOTU)."\n";
		}
	}	
close annFile;

# filte aligned fasta
my $keep = 0;
open(outFasta, ">filte_aligned.fasta");
open fastaFile, $align_fasta;
	while(<fastaFile>){
		#$keep = defined($groupOTU{$1})&&defined($topGenera{$groupOTU{$1}})? 1:0 if($_ =~ /^>(\d+)/);
		$keep = defined($groupOTU{$1})&&defined($topGenera{$groupOTU{$1}})? 1:0 if($_ =~ /^>(.*)/);
		print outFasta $_ if($keep);
	}
close fastaFile;
close outFasta;

# remove '.' character in filte_aligned.fasta
`/app_tools/mothur/mothur "#filter.seqs(processors=10,fasta=filte_aligned.fasta, trump=., vertical=T)"`;

# make phylogenetic tree
`make_phylogeny.py -i filte_aligned.filter.fasta -o phylogenetic_tree.tre`;

# plot
`graphlan_annotate.py --annot annotate_tree.txt phylogenetic_tree.tre phylogenetic_tree.xml`;
`graphlan.py phylogenetic_tree.xml $output_tree --size 14`;
