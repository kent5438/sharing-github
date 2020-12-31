#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);

my $mapping_file= defined $ARGV[0]? shift @ARGV : '';
my $merged_fasta= defined $ARGV[0]? shift @ARGV : '';
my $cluster_list= defined $ARGV[0]? shift @ARGV : '';
my $classify_otu= defined $ARGV[0]? shift @ARGV : '';
my $output_file	= defined $ARGV[0]? shift @ARGV : 'statistic_sample.txt';

die "perl $0 [mapping txt] [merged fasta] [cluster list] [classify otu] ([output file: statistic_sample.txt])\n" if(!-f $mapping_file || !-f $merged_fasta || !-f $cluster_list || !-f $classify_otu);
make_path(dirname($output_file));

# read mapping file
my $idIndex = 0;
my @sample = ();
open mapFile, $mapping_file;
	while(<mapFile>){
		chomp;
		my @column = split(/\t/, $_);
		next if($#column < 3);
		if($_ =~ /^#/){
			for my $i(0..$#column){
				$idIndex = $i if($column[$i] eq '#SampleID');
			}
		}else{
			push @sample, $column[$idIndex];
		}
	}
close mapFile;
@sample = sort @sample;

### statistic tags
my %fasta = %{&ReadFasta($merged_fasta)}; # read merge fasta
my %totalTags = ();
my %uniqTags = ();
for my $i( 0 .. $#sample){
	my @cache = map { m/^$sample[$i]\_\d+.+/g } keys %fasta;
	$totalTags{$sample[$i]} = $#cache +1;
	
	my %uniqSeq = ();
	foreach (@cache){
		$uniqSeq{$fasta{$_}} = 1;
	}
	$uniqTags{$sample[$i]} = keys %uniqSeq;
}



### statistic otus

# read cluster list to group sequence
my @cloumnTitle = ();
my %clusterGroup = ();
my $title = 1;
open(clusterFile, $cluster_list);
    while(<clusterFile>){
        chomp;
        if($title){
            $title = 0;
            @cloumnTitle = split(/\t/, $_);
        }else{
            my @column = split(/\t/, $_);
            for(my $i=0; $i<=$#cloumnTitle; $i++){
                next if($cloumnTitle[$i] eq 'label' || $cloumnTitle[$i] eq 'numOtus');
                $clusterGroup{$cloumnTitle[$i]} = $column[$i];
            }
            last;
        }
    }
close clusterFile;

# read classify otu
my %taxaTags = ();
my %OTU = ();
$title = 1;
open(classifyFile, $classify_otu);
    while(<classifyFile>){
        chomp;
        if($title){
            $title = 0;
            next;
        }
        my ($OTU, $size, $taxonomy) = split(/\t/, $_);
        my @names = split(/,/, $clusterGroup{$OTU});
        
        foreach(@names){
            if($_ =~ /^(.+)\_\d+$/){
                $OTU{$1}++;
                $taxaTags{$1}++ if($taxonomy !~ /^unknown/);
            }
        }
    }
close classifyFile;



### output
open(outFile, ">$output_file");
	print outFile join("\t", 'SampleID', 'Total_Tags', 'Unique_Tags', 'Taxon_Tags', 'OTUs') ."\n"; # print title
	
	for my $i(0.. $#sample){
		# init for none case
		my $name = $sample[$i];
		$totalTags{$name} = 0 if(!defined($totalTags{$name}));
		$uniqTags{$name} = 0 if(!defined($uniqTags{$name}));
		$taxaTags{$name} = 0 if(!defined($taxaTags{$name}));
		$OTU{$name} = 0 if(!defined($OTU{$name}));
		
		print outFile join("\t", $name, $totalTags{$name}, $uniqTags{$name}, $taxaTags{$name}, $OTU{$name}) ."\n";
	}
	
close outFile;

sub ReadFasta(){
	my $fastaFile = shift @_;
	my %sequence;
	my $name = '';
	open(fastaF, $fastaFile);
		while(<fastaF>){
			chomp;
			next if($_ eq '');
			if($_ =~ /^>(\S+)/){
				$name = $1;
			}else{
				$sequence{$name} .= $_;
			}
		}
	close(fastaF);
	return \%sequence;
}