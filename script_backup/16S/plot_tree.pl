#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);

my $tree_file   = defined $ARGV[0]? $ARGV[0] : '';
my $mapping_file= defined $ARGV[1]? $ARGV[1] : '';
my $output_file	= defined $ARGV[2]? $ARGV[2] : 'UPGMA_tree.png';
my $group_title	= defined $ARGV[3]? $ARGV[3] : 'Description';
my $output_title= defined $ARGV[4]? $ARGV[4] : 'UPGMA tree';

die "perl $0 [tree file] [mapping txt] ([output: UPGMA_tree.png]) ([group title: Description]) ([figure title:UPGMA tree])\n" if(!-f $tree_file || !-f $mapping_file);
make_path(dirname($output_file));

# setup ARGV
my $outWidth  = 1080;
my $outHeight = 1080;

# read mapping file for group
my %group = ();
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
		}
	}
close mapFile;

# read tree
my $tree = '';
open treeFile, $tree_file;
	while(<treeFile>){
		chomp;
		$tree .= $_;
	}
close treeFile;

# replace tip label & output cache.tre
my %code = ();
my $maxLabel = 0;
foreach my $name(keys %group){
	my $label = $group{$name};
### modified at 20181112
#	if($group{$name} ne $name){
		$tree =~ s/\b${name}\b/$name\_|_$group{$name}/;
		$label = "$name | $group{$name}";
#	}
###
	$code{$group{$name}} .= "'$label', ";
	$maxLabel = length($label) if(length($label)>$maxLabel);
}

# output cache tree file
open(cacheTree, ">cache.tre");
	print cacheTree $tree;
close cacheTree;

# setup R code
open(cacheR, ">tree.R");
	my $tipMultiple = 1 + $maxLabel*0.02; # Magnification 0.02 per char
	# setup group code
	my $groupCode = 'list(';
	foreach my $groupName (sort {$a cmp $b} keys %code){
		$groupCode .= "'$groupName'=c(". substr($code{$groupName}, 0, -2) ."), ";
	}
	$groupCode = substr($groupCode, 0, -2) .")";
	
	# output R code
	print cacheR << "START";
library('ggtree')
tree <- read.tree('cache.tre')
tree\$tip.label <- gsub("_\\\\|_", " \\\\| ", tree\$tip.label)
groupInfo <- $groupCode
groupTree <- groupOTU(tree, groupInfo)

png('$output_file', width=$outWidth, height=$outHeight)
p <- ggtree(groupTree, aes(color=group), cex=2)
p <- p + ggtitle("$output_title group by $group_title")
p <- p + ggplot2::xlim(0, max(tree\$edge.length) * 2.5)
p <- p + geom_tiplab(size=9) 
p <- p + theme_bw(base_size=30)
p <- p + theme(legend.position="none", plot.title=element_text(hjust=0.5))
print(p)
dev.off()
START

close cacheR;

`Rscript tree.R`;
