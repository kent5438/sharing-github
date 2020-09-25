#!/usr/bin/perl -w
# Author: Daniel Liu
# perl report-ITS.pl.pl --sample=list.txt
# Date: 2016.11.25
# We use HISAT and StringTie to build new gtf file (gene anotation); We don't use Ballgown because it need more replicate in each condition
# We use RSEM for differential Gene Expression, GO/KEGG enrichment
# Data structure
#  
# Function:
#
#

$main = new report;
$main-> report;

package report;
use strict;
use warnings;
use Getopt::Long;
use lib '/export/EC1680U/perl/lib'; # Global
use Daniel::Config qw(init_cfg);
use Daniel::Utility qw(process_cmds_serial process_cmds_parallel process_cmd create_fullshell checkjob checkjobid);
use Daniel::Report qw (createBarcodeHtml createQCHtml createOtuHtml createTaxonomyHtml);
sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;
  $self->init;
  return $self;
}
sub init {
  our %cfg = init_cfg ("report-ITS.cfg");
  our $dir = `pwd`; chomp($dir); # workdir of script
  our $output = "$dir/report";
  my $self = shift;
  ########### template ###########
  our $barcode_template = $cfg{template}{Barcode};
  our $index_template = $cfg{template}{Index};
  our $otustat_template  = $cfg{template}{OtuStat};
  our $readstat_template = $cfg{template}{ReadStat};
  our $taxonomystat_template = $cfg{template}{TaxonomyStat};
  ########### template ###########
  
  ########### file ###########
  our $barcodelist = $cfg{file}{BarcodeList};  
  our $otustat = $cfg{file}{OtuStat};
  our $readstat = $cfg{file}{ReadStat};  
  our $taxonomystat = $cfg{file}{TaxonomyStat};
  ########### file ###########
  
  ########### file ###########
  our $samplesize = $cfg{parameter}{SampleSize};  
  
  ########### file ###########
  
  unless (-e $output) {system ("mkdir -p $output");}
}
sub report {
  our $barcodelist;
  our $barcode_template;
  our $index_template;
  our $otustat_template;
  our $otustat;
  our $output;
  our $readstat;
  our $readstat_template;
  our $samplesize;
  our $taxonomystat;
  our $taxonomystat_template;
  my $self = shift;
  my $content = "";
  my $cmd;
    
  $self->inithtml($output);
  # barcodes statistics
  createBarcodeHtml ($barcodelist, $barcode_template,  "$output/html");
  $content = $content."<li><a href=html/barcode.html target=display>Barcode Stats</a></li>";
  # QC statistics
  createQCHtml ($readstat, $readstat_template,  "$output/html");
  $content = $content."<li><a href=html/readstat.html target=display>Read Stats</a></li>";
  # otu statistics
  createOtuHtml ($otustat, $otustat_template,  "$output/html");
  $content = $content."<li><a href=html/otustat.html target=display>OTU Stats</a></li>";
  # taxonomy statistics
  createTaxonomyHtml ($taxonomystat, $taxonomystat_template,  "$output/html");
  $content = $content."<li><a href=html/taxonomy.html target=display>Taxonomy Stats</a></li>";
  
  #rarefaction_plots
  $content = $content."<li><a href=rank_abundance_graph/rank_abundance_graph.pdf target=display>rank abundance graph</a></li>";
  $content = $content."<li><a href=rare_alpha_plot/rarefaction_plots.html target=display>rarefaction plots</a></li>";
  
  $content = $content."<a style=color: blue; text-align: center; display: block; background-color: transparent; text-decoration: none; href=AlphaBetaDiversity/alphaDiversity.xls target=display>Alpha diversity (.xls)</a><br>";
  $content = $content."<a style=color: blue; text-align: center; display: block; background-color: transparent; text-decoration: none; href=AlphaBetaDiversity/betaDiversity.xls target=display>Beta diversity (.xls)</a>";
  $content = $content."<HR>";
    
  #taxa_plots bar_charts
  $content = $content."<li><a href=taxa_plots/bar_charts.html target=display>taxa_plots bar plots</a></li>";
  
  if ($samplesize >= 3) {
    #taxa_plots area_chartt
    $content = $content."<li><a href=taxa_plots/area_charts.html target=display>taxa_plots area plots</a></li>";
	#OTUs_heatmap
    $content = $content."<li><a href=OTUs_heatmap/OTUs_heatmap.pdf target=display>OTU heatmap (genus)</a></li>";
	#PCA
	$content = $content."<li><a href=files/PCA.pdf target=display>PCA plot</a></li>";
  }
  
  $cmd = "sed 's|LeftList|$content|' $index_template > $output/index.html";
  process_cmd ($cmd);   
}

#system("Rscript $dir/PCA_plot.R $output");

sub inithtml {
  my $self = shift;
  my $output = shift;
  my $cmd;
  my @cmds;
  unless (-e "$output/images") {system("mkdir -p $output/images");} # image files
  unless (-e "$output/html")   {system("mkdir -p $output/html");} # html files	
  unless (-e "$output/files")  {system("mkdir $output/files");}   # files
  $cmd = "cp -r /export/EC1680U/kentchen/Template/ITS/css $output";
  push (@cmds, $cmd);
  $cmd = "cp -r /export/EC1680U/kentchen/Template/ITS/images/* $output/images";
  push (@cmds, $cmd);
  $cmd = "cp -r /export/EC1680U/kentchen/Template/ITS/DataAnalysisWorkflow.html $output/html/";
  push (@cmds, $cmd);
  $cmd = "cp -r /export/EC1680U/kentchen/Template/ITS/help.pdf $output/files/";
  push (@cmds, $cmd);
#  $cmd = "cp -r `ls metaITS_OUTPUT/* | grep -v PCA.pdf` $output";
#  push (@cmds, $cmd);
  $cmd = "rsync -avq metaITS_OUTPUT/* $output";
  push (@cmds, $cmd);
  process_cmds_parallel(@cmds);
}

1;

# Convert PCoA default eps to pdf
#my @all_eps;
#my @eps_gz = `find report/PCoA -name '*eps.gz'`; chomp(@eps_gz);
#foreach my $eps (@eps_gz){
#    my $prefix = `ls $eps | cut -d '.' -f1`; chomp($prefix);
#    if(! -e "$prefix.eps"){system("zcat $prefix.eps.gz > $prefix.eps");}
#    push(@all_eps, "$prefix.eps");
#}
#my $eps_pdf_count = `find  report/PCoA -name '*eps.pdf' | wc -l`; chomp($eps_pdf_count);
#if($eps_pdf_count == 0){system("/export/EC1680U/software/EPStools/eps2pdf @all_eps");}
