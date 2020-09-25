#!/usr/bin/env perl
#
# Copyright belong to Chen Yuan Liu 
# Author: Chen-Yuan Liu
# Date: 09-21, 2016
# Version 2 Using perl module
# The Core for this pipeline is RSEM
#  
# Function:
#
$main = new RNASEQCoreV3;
$main-> parse_command_line;
print 'CalculateExpression start'.$/;
$main-> CalculateExpression;
print 'GenerateDataMatrix start'.$/;
$main-> GenerateDataMatrix;
print 'ebseq start'.$/;
$main-> ebseq;
print 'report start'.$/;
$main-> report;

package RNASEQCoreV3;
use strict;
use warnings;
use lib '/export/EC1680U/perl/lib'; # Horse
use Getopt::Long;
use vars qw (%cfg $parameters @samples);
use Daniel::Config qw(init_cfg);
use Daniel::Utility qw(process_cmds_serial process_cmds_parallel process_cmd create_fullshell checkjob checkjobid);
use Daniel::Report qw (fastqstats procFastqStats createFastqStatsHtml createEnrichHtml procRSEMAlignStats createRSEMAlignStatsHtml );


sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;  
  return $self;
}
sub init {
  our %cfg = init_cfg ("dge.cfg");  
  our $dir = `pwd -L`; chomp($dir); # workdir of script 
  our $output     = "output";
  our $out_CE     = "$dir/$output/1.CalculateExpression";
  our $out_DM     = "$dir/$output/2.GenerateDataMatrix";  
  our $out_ebseq  = "$dir/$output/3.ebseq";
  our $stats_out  = "$dir/$output/4.stats_out";  
  our $report_out = "$dir/$output/5.report";      
  
  my $self = shift;
  ########### Program ########### 
  our $bowtie2_root = $cfg{software}{BOWTIE2_ROOT};
  our $perlbin_root = $cfg{software}{Perlbin_ROOT};
  our $rsem_root = $cfg{software}{RSEM_ROOT};
  our $fastqstats_root= $cfg{software}{FastqStats_ROOT};
  ########### Program ########### 
  
  ########### parameter ########### 
  our $genome = $cfg{parameter}{genome};  
  ########### parameter ########### 
  
  ########### File ########### 
  our $reference = $cfg{files}{REFERENCE};
  our $ngvec     = $cfg{files}{NGVEC};
  our $rscript = $cfg{files}{Rscript};
  ########### File ########### 
  
  ########### template ###########
  our $fastqstats_template = $cfg{template}{FASTQSTATS};
  our $samstats_template = $cfg{template}{SAMSTATS};
  our $enrich_template = $cfg{template}{ENRICH};
  our $index_template = $cfg{template}{INDEX};
  ########### template ###########
  
  unless (-e $out_CE)    {system ("mkdir -p $out_CE");}
  unless (-e $out_DM)    {system ("mkdir -p $out_DM");}
  unless (-e $out_ebseq) {system ("mkdir -p $out_ebseq");}
  unless (-e $stats_out)    {system ("mkdir -p $stats_out");}
  unless (-e $report_out) {system ("mkdir -p $report_out");}
}
sub parse_command_line {
  our @comparsions;
  our $parameters;
  our @samples;
  my $self = shift;
  my @comma;  
  my $ln;
  my $sample;
  my @semicolon;  
  my %unique = ();
  if(scalar(@ARGV) < 1){
    print usage();
    exit(1);
  } 

  GetOptions(
    'sample=s'  => \$sample,
  );
#SET:<case1[,case2,...]>;<control1[,control2,...]>
#NOTE: the label should be exactly the readname 
#SET1:T24-1ug;T24
#SET2:T24-200ng;T24
#SET3:T24-2ug;T24
#SET4:T24-500ng;T24
 
  open (IN, '<', $sample) or die "Can't read $sample file";
  while (<IN>) {
    $ln = $_;
    chomp($ln);
    if ($ln=~/(^SET.*):(.*)/) {
      @semicolon = split(/\;/, $2);	  
	  $parameters->{$1}->{case} = $semicolon[0];
	  @comma = split(/,/, $semicolon[0]);
	  push (@samples, @comma);	  
	  $parameters->{$1}->{ctrl} = $semicolon[1];	
	  @comma = split(/,/, $semicolon[1]);
	  push (@samples, @comma);
	  push (@comparsions, "$semicolon[0]-$semicolon[1]");
    }
  }
  close IN;
  # Get Unique Sample Name
  foreach my $item (@samples) {
    $unique{$item} ++;
  }
  @samples = keys %unique;  
  $self->init;
}
sub usage {
  my $self = shift;
  my $msg = <<_EOUSAGE_;
# CMD:
# perl $0 --sample=list.txt 

[software]
BOWTIE2_ROOT=/export/EC1680U/software/bowtie2-2.2.6
FastqStats_ROOT=/export/EC1680U/software/ExpressionAnalysis-ea-utils-27a4809/clipper
Perlbin_ROOT=/export/EC1680U/perl/bin
RSEM_ROOT=/export/EC1680U/software/RSEM-latest

[files]
REFERENCE=/export/EC1680U/DataBase/hg19/RSEM/reference/hg19.refseq_125polyA
NGVEC=/export/EC1680U/DataBase/hg19/RSEM/reference/hg19.refseq_125polyA.transcripts.ngvec
Rscript=/export/EC1680U/Rscript/dge

[parameter]
genome=human

[template]
FASTQSTATS=/export/EC1680U/kentchen/Template/RnaSeq/fastqstats_template
SAMSTATS=/export/EC1680U/kentchen/Template/RnaSeq/samstats_template
INDEX=/export/EC1680U/kentchen/Template/RnaSeq/index_template
ENRICH=/export/EC1680U/kentchen/Template/RnaSeq/enrich_template
_EOUSAGE_
;
  return $msg;
}
sub CalculateExpression {
  our $bowtie2_root;
  our $dir;
  our @samples;
  our $rsem_root;
  our $rsem_calculate_expression;
  our $reference ;  
  our $out_CE;  
  our $report_out;
  our $perlbin_root;
  my $self = shift;  
  my $cmd;  
  my $err;
  my $log;
  my $qname = "CalculateExpression";
  my $r1;
  my $r2;
  my $sample;
  my $shell;
  foreach $sample (@samples) {    
    $shell = "$out_CE/$sample.CalculateExpression.sh";
	$err   = "$out_CE/$sample.CalculateExpression.err";
	$log   = "$out_CE/$sample.CalculateExpression.log";
    $cmd = "";    
    $r1 = "$dir/reads/$sample.R1.fastq.gz";
#    $r2 = "$dir/reads/$sample.R2.fastq.gz";	
	$cmd=$cmd."export PATH=/export/EC1680U/software/anaconda2/bin:\$PATH\n";
	$cmd=$cmd."export PERL5LIB=/export/EC1680U/software/anaconda2/lib/perl5:\$PERL5LIB\n";
	$cmd=$cmd."export PATH=$bowtie2_root:\$PATH\n";
    $cmd=$cmd."$rsem_root/rsem-calculate-expression --bowtie2 -p 20 --output-genome-bam --phred33-quals $r1 $reference $sample \n";          	
	$cmd=$cmd."mv $sample* $out_CE\n";
#    unless (-e $err) {create_fullshell ($shell, $qname, $err, $log, $cmd);  `qsub $shell >> ./jobs` ;}
    unless (-e $err) {create_fullshell ($shell, $qname, $err, $log, $cmd);  `sh $shell 2>$err 1>$log` ;}
  }
#  checkjobid("");
  
    # remove params (--paired-end --strand-specific --forward-prob 0) to fit BGI's data
    # --forward-prob <double>
        # Probability of generating a read from the forward strand of a
        # transcript. Set to 1 for a strand-specific protocol where all
        # (upstream) reads are derived from the forward strand, 0 for a
        # strand-specific protocol where all (upstream) read are derived from
        # the reverse strand, or 0.3 for a non-strand-specific protocol.
        # (Default: 0.3)
    # --strand-specific
        # The RNA-Seq protocol used to generate the reads is strand specific,
        # i.e., all (upstream) reads are derived from the forward strand. This
        # option is equivalent to --forward-prob=1.0. With this option set, if
        # RSEM runs the Bowtie/Bowtie 2 aligner, the '--norc' Bowtie/Bowtie 2
        # option will be used, which disables alignment to the reverse strand
        # of transcripts. (Default: off)
	
}

sub GenerateDataMatrix {
  our $parameters;
  our $out_CE;
  our $out_DM;  
  our $rsem_root;
  my $self = shift;
  my $case;
  my $case_gen = "";
  my $case_iso = "";
  my @cases;
  my $ctrl;
  my $ctrl_gen = "";
  my $ctrl_iso = "";
  my @ctrls;
  my $cmd;  
  my $err;
  my $idx;  # index of case & control
  my $log;
  my $qname = "GenerateDataMatrix";
  my $set;
  my $shell;
  
# SET1:T24-1ug;T24
  foreach $set (sort keys %{$parameters}) {
    $case = $parameters->{$set}->{case};
	$ctrl = $parameters->{$set}->{ctrl};
	unless (-e "$out_DM/$case-$ctrl") {system ("mkdir -p $out_DM/$case-$ctrl");}
	$shell = "$out_DM/$case-$ctrl/GenerateDataMatrix.sh";
	$err   = "$out_DM/$case-$ctrl/GenerateDataMatrix.err";
	$log   = "$out_DM/$case-$ctrl/GenerateDataMatrix.log";
    $cmd = "";    
	
	@cases = split(/\,/, $case);
	@ctrls = split(/\,/, $ctrl);
	$case_gen = ""; $case_iso = "";
	for ($idx = 0;$idx <= $#cases; $idx++) {
	  $case_gen = $case_gen." $out_CE/$cases[$idx].genes.results ";
	  $case_iso = $case_iso." $out_CE/$cases[$idx].isoforms.results ";
	}
	$ctrl_gen = ""; $ctrl_iso = "";
	for ($idx = 0;$idx <= $#ctrls; $idx++) {
	  $ctrl_gen = $ctrl_gen." $out_CE/$ctrls[$idx].genes.results ";
	  $ctrl_iso = $ctrl_iso." $out_CE/$ctrls[$idx].isoforms.results ";
	}
	$cmd=$cmd."$rsem_root/rsem-generate-data-matrix $case_gen $ctrl_gen | sed 's|$out_CE/||g' > $out_DM/$case-$ctrl/$case-$ctrl.genMat\n";
	$cmd=$cmd."$rsem_root/rsem-generate-data-matrix $case_iso $ctrl_iso | sed 's|$out_CE/||g' > $out_DM/$case-$ctrl/$case-$ctrl.isoMat\n";
#	unless (-e $err) {create_fullshell ($shell, $qname, $err, $log, $cmd);   `qsub $shell >> ./jobs`  ;}
	if(! -e "$out_DM/$case-$ctrl/$case-$ctrl.genMat"){create_fullshell ($shell, $qname, $err, $log, $cmd);   `sh $shell 2>$err 1>$log`  ;}
  }   
# checkjobid("");
}
sub ebseq {  
  our $out_DM;
  our $out_ebseq;
  our $parameters;
  our $ngvec;
  our $rsem_root;
  
  my $self = shift;
  my $case;
  my $caseno;
  my @cases;
  my $cmd;
  my $ctrl;
  my @ctrls;
  my $ctrno;  
  my $err;
  my $log;
  my $qname = "ebseq";
  my $set;
  my $shell;
  
# SET1:T24-1ug;T24
  foreach $set (sort keys %{$parameters}) {
    $case = $parameters->{$set}->{case};
	$ctrl = $parameters->{$set}->{ctrl};
	@cases = split(/\,/, $case);
	@ctrls = split(/\,/, $ctrl);
	$caseno = scalar @cases;
	$ctrno  = scalar @ctrls;   
	unless (-e "$out_ebseq/$case-$ctrl") {system ("mkdir -p $out_ebseq/$case-$ctrl");}
	$shell = "$out_ebseq/$case-$ctrl/ebseq.sh";
	$err   = "$out_ebseq/$case-$ctrl/ebseq.err";
	$log   = "$out_ebseq/$case-$ctrl/ebseq.log";
    $cmd = "";    
	$cmd=$cmd."$rsem_root/rsem-run-ebseq --ngvector $ngvec $out_DM/$case-$ctrl/$case-$ctrl.genMat $caseno,$ctrno $out_ebseq/$case-$ctrl/$case-$ctrl.genMat.result \n";	
	$cmd=$cmd."$rsem_root/rsem-run-ebseq --ngvector $ngvec $out_DM/$case-$ctrl/$case-$ctrl.isoMat $caseno,$ctrno $out_ebseq/$case-$ctrl/$case-$ctrl.isoMat.result \n";	
    $cmd=$cmd."$rsem_root/rsem-control-fdr $out_ebseq/$case-$ctrl/$case-$ctrl.genMat.result 0.05 $out_ebseq/$case-$ctrl/$case-$ctrl.genMat.fdr-0.05.result \n";
	$cmd=$cmd."$rsem_root/rsem-control-fdr $out_ebseq/$case-$ctrl/$case-$ctrl.isoMat.result 0.05 $out_ebseq/$case-$ctrl/$case-$ctrl.isoMat.fdr-0.05.result \n";
#	$cmd=$cmd."$rsem_root/rsem-control-fdr $out_ebseq/$case-$ctrl/$case-$ctrl.genMat.result 0.3 $out_ebseq/$case-$ctrl/$case-$ctrl.genMat.fdr-0.3.result \n";
#	$cmd=$cmd."$rsem_root/rsem-control-fdr $out_ebseq/$case-$ctrl/$case-$ctrl.isoMat.result 0.3 $out_ebseq/$case-$ctrl/$case-$ctrl.isoMat.fdr-0.3.result \n";
#	unless (-e $err) {create_fullshell ($shell, $qname, $err, $log, $cmd);   `qsub $shell >> ./jobs` ;}
	if(! -e "$out_ebseq/$case-$ctrl/$case-$ctrl.genMat.fdr-0.05.result"){create_fullshell ($shell, $qname, $err, $log, $cmd);   `sh $shell 2>$err 1>$log`  ;}
  }     
#  checkjobid("");
}

sub report {
  our @comparsions;
#  my @comparsions;
  our $dir;
  our $enrich_template;
  our $fastqstats_root;
  our $fastqstats_template;  
  our $genome;  
  our $index_template;
  our $out_CE; 
  our $out_DM ;
  our $out_ebseq;
  our $parameters;
  our $perlbin_root;
  our $report_out;   
  our $rscript ;
  our $samstats_template;  
  our $stats_out;  
  our @samples;
  my $self = shift;  
  my @FPKMs;
  my $case; 
  my @comma;
  my $cmd;
  my $comparsion;
  my $content;
  my $ctrl;    
  my $data; 
  my $err;
  my $genMat;
  my $log;
  my $output; 
  my $qname = "report";
  my $r1;
  my $r2;
  my $sample;
  my $set;
  my $shell;
  my $tmp; # temp string
 
  $self->inithtml ($report_out);

# All samples heatmap plot by merged FPKM	
	foreach $sample (@samples){
        open IN, "<$out_CE/$sample.genes.results";
        open OUT, ">$out_CE/$sample.FPKM.txt";
        while(<IN>){
            chomp;
            $_ =~ s/FPKM/$sample/g;
            my @sp = split/\t/, $_;
            print OUT "$sp[0]\t$sp[6]\n";
        }
        close IN;
        close OUT;
        push(@FPKMs, "$out_CE/$sample.FPKM.txt");
    }

    my $merge_cmd;
    my @list;
    my $num = 1;
    foreach my $FPKM (@FPKMs){
        my $name = "tmp$num";
        $merge_cmd .= $name.' <- read.table("'.$FPKM.'", header=T, check.names=F, fill=TRUE, quote="\"", sep="\t");';
        push(@list, $name);
        $num++
    }
    my $list_comma = join(',', @list);
    $merge_cmd .= 'merge <- Reduce(function(x,y) merge(x, y, by="gene_id", all=TRUE), list('.$list_comma.'));';
    $merge_cmd .= 'merge <- merge[, !duplicated(colnames(merge))];';
    $merge_cmd .= 'write.table(merge, file="'.$out_CE.'/all_samples_FPKM.txt", sep="\t", quote=F, na="-", row.names = FALSE);';
    system("Rscript" . " -e '" . $merge_cmd . "'");

    system("Rscript /export/EC1680U/Rscript/dge/heatmap_eachSample.R $out_CE/all_samples_FPKM.txt $report_out/images");
    system("$perlbin_root/tab2xls.pl $out_CE/all_samples_FPKM.txt $report_out/files/all_samples_FPKM.xls");
	
# All samples correlation by merged FPKM
	if (! -e "$out_CE/all_samples_FPKM.txt"){die "### ERROR: Please make sure $out_CE/all_samples_FPKM.txt is generated!\n";}
	system("Rscript /export/EC1680U/Rscript/dge/correlation_allSample.R $out_CE/all_samples_FPKM.txt $report_out/images/Sample_correlation.pdf");


  foreach $set (sort keys %{$parameters}) {
    $cmd = "";
    @samples = ();
    $case = $parameters->{$set}->{case};
    $ctrl = $parameters->{$set}->{ctrl};	    
    $output = "$report_out/$case-$ctrl";
    unless (-e "$output") {system ("mkdir -p $output");}  
    
    $shell = "$output/report.sh";
    $log   = "$output/report.log";
    $err   = "$output/report.err";
    @comma = split(/,/, $case);	
    push (@samples, @comma);
    @comma = split(/,/, $ctrl);
    push (@samples, @comma);
#    push (@comparsions, "$case-$ctrl");
    
    $genMat = "$out_DM/$case\-$ctrl/$case-$ctrl.genMat";
      
    # expression level of genes  [read count]
    $cmd=$cmd."awk -F \" \" '{if (\$2 ~/results/) {print \"Gene\" \$0} else {print \$0}}' $genMat | sed 's|.genes\\.results||g' | sed 's|\"||g'> $output/$case-$ctrl.genMat \n" ; 	  	     
    $cmd=$cmd."$perlbin_root/tab2xls.pl $output/$case-$ctrl.genMat $output/$case-$ctrl.genMat.xls\n";

    # FPKM 
    $cmd=$cmd."$perlbin_root/FPKM.pl --cases=$case --ctrls=$ctrl --input=$out_CE --output=$output \n";	
    
    # Differential expression of RSEM transcript
    $cmd=$cmd."sed 's|\"||g' $out_ebseq/$case-$ctrl/$case-$ctrl.genMat.result | sed 's|PPEE\\tPPDE\\tPostFC\\tRealFC\\tC1Mean\\tC2Mean|Gene\\tPPEE\\tPPDE\\tPostFC\\tRealFC\\tC1Mean\\tC2Mean|' > $out_ebseq/$case-$ctrl/$case-$ctrl.genMat.result.1 \n";	
    $cmd=$cmd."$perlbin_root/tab2xls.pl $out_ebseq/$case-$ctrl/$case-$ctrl.genMat.result.1 $output/$case-$ctrl.genMat.result.xls \n";
    	
    # HierarchicalClustering of read count [gene] 
    # The rscript may fail if there are only two samples
#    if (scalar @samples > 2) {
#      $cmd=$cmd."Rscript $rscript/HierarchicalClustering.r $out_DM/$case-$ctrl/$case-$ctrl.genMat $output\n";	  
#    }
    
    # Scatter plot of RSEM transcript expression
    $cmd=$cmd."Rscript $rscript/Scatter.new.r $out_ebseq/$case-$ctrl/$case-$ctrl.genMat.result.1 $output \n";
    
      
    # enrichment analysis GO/KEGG for human, mouse, rat
	my $check_GO = `ls $output | grep "enrichGO" | wc -l`; chomp($check_GO);
	my $check_KEGG = "$output/keggmap";
	if($check_GO == 0 || (! -e $check_KEGG)) {
	    $cmd=$cmd."/usr/bin/Rscript $rscript/clusterProfiler.over-representation.v3.1.kentchen.r $out_ebseq/$case-$ctrl/$case-$ctrl.genMat.fdr-0.05.result $genome $output \n";
	#	$cmd=$cmd."Rscript $rscript/clusterProfiler.over-representation.v3.1.kentchen.r $out_ebseq/$case-$ctrl/$case-$ctrl.genMat.fdr-0.3.result $genome $output \n";
	}

    $cmd=$cmd."$perlbin_root/tab2xls.pl $output/gen.enrichKEGG.txt $output/gen.enrichKEGG.xls \n";	
    $cmd=$cmd."$perlbin_root/tab2xls.pl $output/gen.bp.enrichGO.txt $output/gen.bp.enrichGO.xls \n";	
    $cmd=$cmd."$perlbin_root/tab2xls.pl $output/gen.mf.enrichGO.txt $output/gen.mf.enrichGO.xls \n";	
    $cmd=$cmd."$perlbin_root/tab2xls.pl $output/gen.cc.enrichGO.txt $output/gen.cc.enrichGO.xls \n";	
    unless (-e "$output/enrich.html") {createEnrichHtml ($enrich_template, "$output", "./");}	

    $tmp = "";
    foreach $sample (@samples) {
      $tmp=$tmp.$sample."#";
    }	  	
    
    # fastq statistic ; run fastq-stats of ea-utils to calc fq statistic use SGE
    $cmd=$cmd."perl $perlbin_root/fastq-stats.pl --samplelist=$tmp --readdir=$dir/reads --fastqstats_root=$fastqstats_root --stats_out=$stats_out --fastqstats_template=$fastqstats_template --outputdir=$output\n";		
    # process the output of rsem-calculate-expression alignment (bowtie2) !! Need to pass array reference to a subroutine of a module
    $cmd=$cmd."perl $perlbin_root/rsem-stats.pl --samplelist=$tmp  --calexpdir=$out_CE --samstats_template=$samstats_template --single-end --outputdir=$output\n";	

#    unless (-e $err) {create_fullshell ($shell, $qname, $err, $log, $cmd);`qsub $shell >> ./jobs` ;}
	my $check_report = `ls $output | wc -l`; chomp($check_report);
	if(($check_report != 26) or ($check_report != 28)) {create_fullshell ($shell, $qname, $err, $log, $cmd);`sh $shell 2>$err 1>$log` ;}
    
    $content = $content."<div class=dropdown>";
    $content = $content."<li><a class=dropbtn> [$case] vs [$ctrl] </a></li>";
    $content = $content."<div class=dropdown-content>";
    $content = $content."<li><a href=$case-$ctrl/fastqstats.html target=display>Fastq Stats</a></li>";
    $content = $content."<li><a href=$case-$ctrl/alignstats.html target=display>Alignment Stats</a></li>";
    $content = $content."<li><a href=$case-$ctrl/$case-$ctrl.genMat.xls target=_blank>Gene Expression</a></li>";
    $content = $content."<li><a href=$case-$ctrl/fpkm.xls target=_blank>FPKM</a></li>";	  		
    $content = $content."<li><a href=$case-$ctrl/$case-$ctrl.genMat.result.xls target=_blank>Differential Gene Expression</a></li>";
    
#    if (scalar @samples > 2) {
#      $content = $content."<li><a href=$case-$ctrl/HC.png target=_blank>Hierarchical Clustering</a></li>";	  
#    }
    $content = $content."<li><a href=$case-$ctrl/Scatter.html target=_blank>Scatter Plot</a></li>";	  
    $content = $content."<li><a href=$case-$ctrl/enrich.html target=_blank>GO/KEGG Plot</a></li>"; 
    $content = $content."</div>";
    $content = $content."</div>";
    }
	
	  # 計算組間的 heatmap 
  # 
  $tmp = "";
  if (scalar @comparsions > 1) {
    foreach $comparsion (@comparsions) {
      $tmp = $tmp."$out_ebseq/$comparsion/$comparsion.genMat.result.1#";
    }    
    chop ($tmp);	
    $cmd = "perl /export/EC1680U/kentchen/DGE/heatmap.kentchen.pl --comps=$tmp --output=$stats_out";
    process_cmd($cmd);
    $cmd = "perl /export/EC1680U/kentchen/DGE/heatmap_union.pl --comps=$tmp --output=$stats_out";
    process_cmd($cmd);
    
	system("sed -i 's/\.genMat\.result\.1//g' $stats_out/comparsion.txt");
	system("sed -i 's/\.genMat\.result\.1//g' $stats_out/comparsion_union.txt");
    		
    $cmd = "Rscript /export/EC1680U/kentchen/DGE/heatmap_rowname.r $stats_out/comparsion.txt $report_out/images/";
    process_cmd($cmd);
    $cmd = "Rscript /export/EC1680U/kentchen/DGE/heatmap_union.r $stats_out/comparsion_union.txt $report_out/images/";
    process_cmd($cmd);
    $cmd = "$perlbin_root/tab2xls.pl $stats_out/comparsion.txt $report_out/files/comparsion.xls";
    process_cmd($cmd);
    $cmd = "$perlbin_root/tab2xls.pl $stats_out/comparsion_union.txt $report_out/files/comparsion_union.xls";
    process_cmd($cmd);

# include samples heatmap & samples correlation into single line html code
	$cmd = "sed 's|Heatmap_Samples|<li><a href=images/all_samples_heatmap.pdf target=display>Samples HeatMap</a><a style=\"color: blue; background-color: transparent; text-decoration: none;\" href=files/all_samples_FPKM.xls target=_blank>(.xls)</a></li><li><a href=images/Sample_correlation.pdf target=display>Sample Correlation</a></li>|' $index_template > $stats_out/index_template\.1";
	process_cmd ($cmd);
    my $check_comp = `cat $stats_out/comparsion.txt | wc -l`; chomp($check_comp);
    if($check_comp != 1){
        $cmd = "sed -i 's|Heatmap_intersection|<li><a href=images/heatmap.pdf target=display>Groups HeatMap (Inter)</a><a style=\"color: blue; background-color: transparent; text-decoration: none;\" href=files/comparsion.xls target=_blank>(.xls)</a></li>|' $stats_out/index_template\.1";
        process_cmd ($cmd);
    } else {
        system("cp /export/EC1680U/kentchen/Template/DGE/heatmap.html $report_out/images");
        $cmd = "sed -i 's|Heatmap_intersection|<li><a href=images/heatmap.html target=display>Groups HeatMap (Inter)</a><a style=\"color: blue; background-color: transparent; text-decoration: none;\" href=files/comparsion.xls target=_blank>(.xls)</a></li>|' $stats_out/index_template\.1";
    }
    $cmd = "sed -i 's|Heatmap_union|<li><a href=images/heatmap_union.pdf target=display>Groups HeatMap (Union)</a><a style=\"color: blue; background-color: transparent; text-decoration: none;\" href=files/comparsion_union.xls target=_blank>(.xls)</a></li>|' $stats_out/index_template\.1";
    process_cmd ($cmd);
  } else {
	$cmd = "sed 's|Heatmap_Samples||' $index_template > $stats_out/index_template\.1";
	print $cmd."\n";
	process_cmd ($cmd);
    $cmd = "sed -i 's|Heatmap_intersection||' $stats_out/index_template\.1";
	print $cmd."\n";
	process_cmd ($cmd);
    $cmd = "sed -i 's|Heatmap_union||' $stats_out/index_template\.1";
	print $cmd."\n";
    process_cmd ($cmd);
  }  
  #exit();
  
#  checkjobid("");
  $cmd = "sed 's|LeftList|$content|' $stats_out/index_template\.1 > $report_out/index.html";			
  process_cmd ($cmd);
}
sub inithtml {
  my $self = shift;
  my $output = shift;
  my $cmd;
  my @cmds;
  unless (-e "$output/images") {system("mkdir -p $output/images");} # image files
  unless (-e "$output/html")   {system("mkdir -p $output/html");} # html files	
  unless (-e "$output/files")  {system("mkdir $output/files");}   # files
  $cmd = "cp -r /export/EC1680U/kentchen/Template/DGE/css $output";
  push (@cmds, $cmd);
  $cmd = "cp -r /export/EC1680U/kentchen/Template/DGE/images/* $output/images";
  push (@cmds, $cmd);
  $cmd = "cp -r /export/EC1680U/kentchen/Template/DGE/workflow.html $output/html";
  push (@cmds, $cmd);
  $cmd = "cp -r /export/EC1680U/kentchen/Template/Help/DGE_Help.pdf $output/html/help.pdf";
  push (@cmds, $cmd);
  $cmd = "cp -r /export/EC1680U/kentchen/Template/DGE/ExperimentWorkflow.html $output/html";
  push (@cmds, $cmd);
  #$cmd = "cp -r /export/EC1680U/kentchen/Template/RnaSeq/ExperimentWorkflow.html $output/html/";
  #push (@cmds, $cmd);
  process_cmds_parallel(@cmds);
}

1;

# Generate multiQC data
if(! -l "QC.pl"){system("ln -s ../QC.pl .");}
system("perl QC.pl");

# Generate Final full table
if(! -l "fullTable.pl"){system("ln -s ../fullTable.pl .");}
system("perl fullTable.pl");
system("mv output/5.report/files/FullTable.xlsx output/5.report");
#system("/export/EC1680U/perl/bin/tab2xls.pl FullTable.txt output/5.report/FullTable.xls");
