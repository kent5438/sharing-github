#!/usr/bin/perl
use strict;
use warnings;
use File::Path qw(make_path);

my  $input_taxa   = defined $ARGV[0]? $ARGV[0] : '';
my  $mapping_file = defined $ARGV[1]? $ARGV[1] : '';
our $level        = defined $ARGV[2]? $ARGV[2] : 'Species';
our $output_folder= defined $ARGV[3]? $ARGV[3] : '3A';

die "perl $0 [input taxa txt] [mapping txt] ([level: Species]) ([output folder: 3A])\n" if(!-f $input_taxa || !-f $mapping_file || $level =~ /^\s?+$/);
make_path($output_folder);

my $qiimeScript = '/usr/local/lib/python2.7/dist-packages/qiime/support_files/R/loaddata.r';

# read mapping file for environment
my %env = ();
my %envIndex = ();
my $idIndex = 0;
open mapFile, $mapping_file;
	while(<mapFile>){
		chomp;
		my @column = split(/\t/, $_);
		next if($#column < 3);
		if($_ =~ /^#/){
			for my $i(0..$#column){
				if($column[$i] eq '#SampleID'){
					$idIndex = $i;
				}elsif($column[$i] =~ /^env_(.+)$/){
					$env{$1} = {};
					$envIndex{$1} = $i;
				}
			}
		}else{
			foreach my $key(keys %envIndex){
				$env{$key}{$column[$idIndex]} = $column[$envIndex{$key}];
			}
		}
	}
close mapFile;

# setup R code
my $code = '';
open(Rcode, ">$output_folder/plot_3A.R");
	# setup environment code
	my $envCode = '';
	my $dataFrame = '';
	if(keys %env){
		foreach my $name (sort {$a cmp $b} keys %env){
			my %envInfo = %{$env{$name}};
			
			my $infoCode = '';
			foreach my $sample(sort {$a cmp $b} keys %envInfo){
				$infoCode .= "'$sample'=$envInfo{$sample}, "
			}
			$infoCode =~ s/,\s$//g;
			
			$envCode .= "$name <- c($infoCode)\n";
			$dataFrame .= "$name, ";
		}
		$dataFrame =~ s/,\s$//g;
		$envCode .= "envData <- data.frame($dataFrame)\n";
	}
		
	# setup basic code
	$code .= "library('vegan')\n";
	$code .= "source('$qiimeScript')\n";
	$code .= "otu.table <- load.qiime.otu.table('$input_taxa')\n";
	if($envCode ne ''){
		$code .= $envCode;
		$envCode = 'envData';
	};
	
	# plot 3A part
	$code .= &plot_DCA('otu.table');
	$code .= &plot_CCA('otu.table', $envCode);
	$code .= &plot_RDA('otu.table', $envCode);

	print Rcode $code;	
close Rcode;

# plot it
`Rscript $output_folder/plot_3A.R`;


sub plot_DCA(){
	my $table = shift @_;
	my $code  = '';
	
	$code .= "dca <- decorana($table)\n";
	
	# output - none label
	$code .= "png('$output_folder/DCA.png', width=10, height=10, units='in', res=300)\n";
	$code .= "par(oma=c(0, 0, 0, 2))\n";
	$code .= "plot(dca, type='p', main='DCA', cex.lab=0.7, cex.axis=0.7, display='sites')\n";
	$code .= "legend('topright', legend='Samples', col=1, pch=1, cex=0.7)\n";
	$code .= "dev.off()\n";	
	
	# output - label
        $code .= "png('$output_folder/DCA_label.png', width=10, height=10, units='in', res=300)\n";
        $code .= "par(oma=c(0, 0, 0, 2))\n";
        $code .= "plot(dca, type='p', main='DCA', cex.lab=0.7, cex.axis=0.7, display='sites')\n";
        $code .= "legend('topright', legend='Samples', col=1, pch=1, cex=0.7)\n";
	$code .= "text(dca, pos=1, cex=0.5)\n";
	$code .= "dev.off()\n";

	return $code;
}

sub plot_CCA(){
	my $table = shift @_;
	my $envData =  shift @_;
	my $code  = '';
	
	$code .= ($envData eq '')? "ccaData <- cca($table)\n" : "ccaData <- cca($table, $envData)\n";
	
	# output - none label
	$code .= "png('$output_folder/CCA.png', width=10, height=10, units='in', res=300)\n";
	$code .= "par(oma=c(0, 0, 0, 2))\n";
	$code .= "plot(ccaData, type='p', main='CCA', xlab='CCA1', ylab='CCA2', cex.lab=0.7, cex.axis=0.7, display=c('sites','wa','bp','cn'))\n";
	$code .= "legend('topright', legend='Samples', col=1, pch=1, cex=0.7)\n";
	$code .= "dev.off()\n";	
	
	# output - label
        $code .= "png('$output_folder/CCA_label.png', width=10, height=10, units='in', res=300)\n";
        $code .= "par(oma=c(0, 0, 0, 2))\n";
        $code .= "plot(ccaData, type='p', main='CCA', xlab='CCA1', ylab='CCA2', cex.lab=0.7, cex.axis=0.7, display=c('sites','wa','bp','cn'))\n";
        $code .= "legend('topright', legend='Samples', col=1, pch=1, cex=0.7)\n";
	$code .= "text(ccaData, pos=1, cex=0.5)\n";
	$code .= "dev.off()\n";
	
	return $code;
}

sub plot_RDA(){
	my $table = shift @_;
	my $envData =  shift @_;
	my $code  = '';
	
	$code .= ($envData eq '')? "rdaData <- rda($table)\n" : "rdaData <- rda($table, $envData)\n";
	
	# output - none label
	$code .= "png('$output_folder/RDA.png', width=10, height=10, units='in', res=300)\n";
	$code .= "par(oma=c(0, 0, 0, 2))\n";
	$code .= "plot(rdaData, type='p', main='RDA', xlab='RDA1', ylab='RDA2', cex.lab=0.7, cex.axis=0.7, display=c('sites','wa','bp','cn'))\n";
	$code .= "legend('topright', legend='Samples', col=1, pch=1, cex=0.7)\n";
	$code .= "dev.off()\n";	
	
	# output - label
        $code .= "png('$output_folder/RDA_label.png', width=10, height=10, units='in', res=300)\n";
        $code .= "par(oma=c(0, 0, 0, 2))\n";
        $code .= "plot(rdaData, type='p', main='RDA', xlab='RDA1', ylab='RDA2', cex.lab=0.7, cex.axis=0.7, display=c('sites','wa','bp','cn'))\n";
        $code .= "legend('topright', legend='Samples', col=1, pch=1, cex=0.7)\n";
	$code .= "text(rdaData, pos=1, cex=0.5)\n";
	$code .= "dev.off()\n";
	
	return $code;
}
