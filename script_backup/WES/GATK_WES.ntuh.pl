#!/usr/bin/perl
#
# Copyright belong to Chen Yuan Liu 
# Author: Chen-Yuan Liu
# Date: 01-13, 2017
# Data structure
# ./GATK_HC.pl  --sample=list.txt
# This script is used for one single sample  (HaplotypeCaller), not `-ERC GVCF` cohort analysis workflow, for 胡務亮專案
# One exome sample can't run VQSR

$main = new GATK_HC;
$main-> parse_command_line;
$main-> alignment;
$main-> indelrealignment;
$main-> bqsr;
$main-> gatk;
#$main-> vqsr;
$main-> vep;
$main-> report;

package GATK_HC;
use strict;
use warnings;
use lib '/export/EC1680U/perl/lib'; # Global
use Getopt::Long;
use threads;
use vars qw (%cfg $sample);
use Daniel::Config qw(init_cfg);
use Daniel::Utility qw(process_cmds_serial process_cmds_parallel process_cmd create_fullshell checkjob checkjobid);
use Daniel::Report qw (fastqstats procFastqStats createFastqStatsHtml samstats procSamStats createSamStatsHtml targetcoverage);
sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;
  return $self;
}
sub init {
  our %cfg = init_cfg ("GATK.cfg");
  our $dir = `pwd`; chomp($dir); # workdir of script
  our $qualdata;                 # quality structure
  our $output                = "output";
  our $alignment             = "$dir/$output/1.alignment";
  our $indelrealignment      = "$dir/$output/2.indelrealignment";
  our $bqsr                  = "$dir/$output/3.bqsr";
  our $hc                    = "$dir/$output/4.hc";
  our $vqsr                  = "$dir/$output/5.vqsr";
  our $vep                   = "$dir/$output/6.vep";  
  our $stats_out             = "$dir/$output/7.stats_out";
  our $report                = "$dir/$output/8.report"; 
  my $self = shift;
    
  ########### Program ########### 
  our $bedtools_root   = $cfg{software}{Bedtools_ROOT};
  our $bwa_root        = $cfg{software}{Bwa_ROOT};
  our $fastqstats_root = $cfg{software}{FastqStats_ROOT};
  our $gatk_root       = $cfg{software}{Gatk_ROOT};
  our $q30             = $cfg{software}{q30};
  our $picard_root     = $cfg{software}{Picard_ROOT};
  our $perlbin_root    = $cfg{software}{Perlbin_ROOT};
  our $samtools_root   = $cfg{software}{Samtools_ROOT};
  our $vep_root        = $cfg{software}{Vep_ROOT};    
  ########### Program ########### 
  
  ########### File ###########  
  our $cosmic     = $cfg{file}{Cosmic}; 
  our $dbsnp      = $cfg{file}{Dbsnp};  
  our $hapmap     = $cfg{file}{Hapmap};
  our $genome     = $cfg{file}{Genome};
  our $interval   = $cfg{file}{Interval};
  our $kg_indel   = $cfg{file}{Kg_indel};
  our $kg_snp     = $cfg{file}{Kg_snp};
  our $mill_indel = $cfg{file}{Mill_indel}; 
  our $omni       = $cfg{file}{Omni};
  our $rscript    = $cfg{file}{Rscript}; 
  our $bed        = $cfg{file}{BED}; 
  ########### File ########### 
  
  ########### parameter ########### 

  ########### parameter ########### 
  
  ########### template ###########
  our $fastqstats_template = $cfg{template}{FASTQSTATS};
  our $index_template = $cfg{template}{INDEX};  
  our $samstats_template = $cfg{template}{SAMSTATS};
  our $qc_template = $cfg{template}{QC};
  ########### template ###########
  
  unless (-e $alignment)        {system ("mkdir -p $alignment");}
  unless (-e $indelrealignment) {system ("mkdir -p $indelrealignment");}
  unless (-e $bqsr)             {system ("mkdir -p $bqsr");}
  unless (-e $hc)               {system ("mkdir -p $hc");}
  unless (-e $vqsr)             {system ("mkdir -p $vqsr");}
  unless (-e $vep)              {system ("mkdir -p $vep");}  
  unless (-e $stats_out)        {system ("mkdir -p $stats_out");}
  unless (-e $report)           {system ("mkdir -p $report");}  
}
sub parse_command_line {
  our $sample;
  my $self = shift;
  my $file;
  my $ln;

  if(scalar(@ARGV) < 1){
    print usage();
    exit(1);
  } 
  GetOptions(
    'sample=s'  => \$file,
  );

  open (IN, '<', $file) or die "Can't read $file file";
  while (<IN>) {
    $ln = $_;
    chomp($ln);
    $sample = $ln;
  }
  close IN;  
  $self->init; 
}
sub usage {
  my $self = shift;
  my $msg = <<_EOUSAGE_;
# CMD:
# perl $0 --sample=list.txt 
# REQUIRED:
GATK.cfg
[software]
Bwa_ROOT=/export/EC1680U/software/bwa-0.7.12
FastqStats_ROOT=/export/EC1680U/software/ExpressionAnalysis-ea-utils-27a4809/clipper
Gatk_ROOT=/export/EC1680U/software/GenomeAnalysisTK-nightly-2016-04-30-gce50d08
Picard_ROOT=/export/EC1680U/software/picard-tools-1.84
Perlbin_ROOT=/export/EC1680U/perl/bin
Samtools_ROOT=/export/EC1680U/software/samtools-1.3
Vep_ROOT=/export/EC1680U/software/ensembl-tools-release-81
Bedtools_ROOT=/export/EC1680U/software/bedtools-2.17.0/bin

[file]
Cosmic=/export/EC1680U/DataBase/hg19/resource/Cosmic.hg19.vcf
Dbsnp=/export/EC1680U/DataBase/hg19/resource/dbsnp_138.hg19.vcf
Genome=/export/EC1680U/DataBase/hg19/BWAIndex/hg19.fa
Interval=/export/EC1680U/DataBase/hg19/mutect/S04380110_Regions.interval.list
Kg_indel=/export/EC1680U/DataBase/hg19/resource/1000G_phase1.indels.hg19.sites.vcf
Mill_indel=/export/EC1680U/DataBase/hg19/resource/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
Rscript=/export/EC1680U/Rscript/WES
BED=/export/EC1680U/Pipeline/WES/resource/S04380110_Regions.bed

[parameter]

[template]
FASTQSTATS=/export/EC1680U/Pipeline/Template/WES/fastqstats_template
INDEX=/export/EC1680U/Pipeline/Template/WES/index_template
SAMSTATS=/export/EC1680U/Pipeline/Template/WES/samstats_template
_EOUSAGE_
;
  return $msg; 
}
sub alignment {
  our $alignment;
  our $bwa_root;
  our $dir;
  our $genome;
  our $parameters;
  our $picard_root;
  our $sample;
  my $self = shift;
  my $cmd;
  my $err;  
  my $group;
  my $log;
  my $qname = "alignment";
  my $r1;
  my $r2;
  my $shell;
    
  $cmd = "";
  $shell = "$alignment/$sample.alignment.sh";
  $err   = "$alignment/$sample.alignment.err";
  $log   = "$alignment/$sample.alignment.log";    	
  $r1 = "$dir/reads/$sample.R1.fastq.gz";
  $r2 = "$dir/reads/$sample.R2.fastq.gz";
  $group = "\@RG\\tID:$sample\\tSM:$sample\\tPL:Illumina";  
  # bwa	
  $cmd=$cmd."$bwa_root/bwa mem -M -R \'$group\' -t 20 $genome $r1 $r2 > $alignment/$sample.sam\n";       
  $cmd=$cmd."cat $alignment/$sample.sam| sed 's/\@PG.*/\@PG\tID:bwa\tPN:bwa/' > $alignment/$sample.sam.1\n";  
  $cmd=$cmd."mv $alignment/$sample.sam.1 $alignment/$sample.sam\n";
  # sortsam
  $cmd=$cmd."java -Xmx4g -jar $picard_root/SortSam.jar INPUT=$alignment/$sample.sam OUTPUT=$alignment/$sample.sorted.bam SORT_ORDER=coordinate TMP_DIR=./tmp \n";
  # markduplicate
  $cmd=$cmd."java -Xmx4g -jar $picard_root/MarkDuplicates.jar INPUT=$alignment/$sample.sorted.bam OUTPUT=$alignment/$sample.sorted.dedup.bam METRICS_FILE=$alignment/$sample-metrics.txt \n";
  # buildbamindex
  $cmd=$cmd."java -Xmx4g -jar $picard_root/BuildBamIndex.jar INPUT=$alignment/$sample.sorted.dedup.bam \n";	
  unless (-e $err) {create_fullshell ($shell, $qname, $err, $log, $cmd); `qsub $shell >> ./jobs` ;}
  
  checkjobid("");
}
sub indelrealignment {
  our $alignment;
  our $gatk_root;
  our $genome;  
  our $indelrealignment;
  our $interval;  
  our $kg_indel;  
  our $mill_indel;  
  our $sample;
  my $self = shift;
  my $cmd;
  my $err;
  my $qname = "indelrealignment";
  my $log;  
  my $r1;
  my $r2;
  my $set;
  my $shell;    
  $cmd = "";
  $shell = "$indelrealignment/$sample.indelrealignment.sh";
  $err   = "$indelrealignment/$sample.indelrealignment.err";
  $log   = "$indelrealignment/$sample.indelrealignment.log";    	
  $cmd=$cmd."java -Xmx4g -jar $gatk_root/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $genome -L $interval -I $alignment/$sample.sorted.dedup.bam -known $kg_indel -known $mill_indel -o $indelrealignment/target_intervals.list -nct 1 -nt 24\n";    
  $cmd=$cmd."java -Xmx4g -jar $gatk_root/GenomeAnalysisTK.jar -T IndelRealigner -R $genome -L $interval -I $alignment/$sample.sorted.dedup.bam -known $kg_indel -known $mill_indel -targetIntervals $indelrealignment/target_intervals.list -o $indelrealignment/$sample.sorted.dedup.realign.bam -nct 1 -nt 1\n";    	  
  unless (-e $err) {create_fullshell ($shell, $qname, $err, $log, $cmd);  `qsub $shell >> ./jobs` ;} 
  checkjobid("");  
}
sub bqsr {
  our $bqsr;
  our $dbsnp;
  our $gatk_root;
  our $genome;  
  our $indelrealignment;
  our $interval;  
  our $kg_indel;  
  our $mill_indel;  
  our $parameters;
  my $self = shift;
  my $cmd;
  my $err;
  my $qname = "bqsr";
  my $log;  
  my $r1;
  my $r2;
  my $shell;

  $cmd = "";
  $shell = "$bqsr/bqsr.sh";
  $err   = "$bqsr/bqsr.err";
  $log   = "$bqsr/bqsr.log";    	
  $cmd=$cmd."java -Xmx4g -jar $gatk_root/GenomeAnalysisTK.jar -T BaseRecalibrator -R $genome -L $interval -I $indelrealignment/$sample.sorted.dedup.realign.bam -knownSites $dbsnp -knownSites $kg_indel -knownSites $mill_indel -o $bqsr/$sample.recal_data.table -nct 8 -nt 1 \n";
  $cmd=$cmd."java -Xmx4g -jar $gatk_root/GenomeAnalysisTK.jar -T PrintReads -R $genome -L $interval -I $indelrealignment/$sample.sorted.dedup.realign.bam -BQSR $bqsr/$sample.recal_data.table -o $bqsr/$sample.sorted.dedup.realign.recal.bam -nct 8 -nt 1 \n";	
  unless (-e $err) {create_fullshell ($shell, $qname, $err, $log, $cmd);  `qsub $shell >> ./jobs` ;}
  checkjobid("");  
}	
sub gatk {
  our $bqsr;
  our $dbsnp;
  our $hc;
  our $gatk_root;
  our $genome;
  our $interval;
  our $sample;
  my $self = shift;	
  my $cmd;
  my $err;
  my $qname = "gatk";
  my $log;  
  my $shell;  
  $cmd = "";
  $shell = "$hc/gatk.sh";
  $err   = "$hc/gatk.err";
  $log   = "$hc/gatk.log";   
  $cmd=$cmd."java -Xmx4g -jar $gatk_root/GenomeAnalysisTK.jar -T HaplotypeCaller -R $genome -L $interval -I $bqsr/$sample.sorted.dedup.realign.recal.bam -bamout $hc/$sample.bamout.bam -nct 1 --dbsnp $dbsnp -o $hc/$sample.vcf \n";
  unless (-e $err) {create_fullshell ($shell, $qname, $err, $log, $cmd);  `qsub $shell >> ./jobs` ;}
  checkjobid("");  
}
sub vqsr {
#https://www.broadinstitute.org/gatk/guide/article?id=2805
  our $dbsnp;
  our $hc;
  our $gatk_root;
  our $genome;
  our $hapmap;
  our $interval;
  our $kg_snp;
  our $omni;
  our $mill_indel;
  our $sample;
  our $vqsr;
  my $self = shift;	
  my $cmd;
  my $err;
  my $qname = "vqsr";
  my $log;  
  my $shell;
  
  $cmd = "";
  $shell = "$vqsr/vqsr.sh";
  $err   = "$vqsr/vqsr.err";
  $log   = "$vqsr/vqsr.log";   
  $cmd=$cmd."java -Xmx4g -jar $gatk_root/GenomeAnalysisTK.jar -T VariantRecalibrator -R $genome -input $hc/$sample.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap -resource:omni,known=false,training=true,truth=true,prior=12.0 $omni -resource:1000G,known=false,training=true,truth=false,prior=10.0 $kg_snp -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $vqsr/recalibrate_SNP.recal -tranchesFile $vqsr/recalibrate_SNP.tranches -rscriptFile $vqsr/recalibrate_SNP_plots.R \n";
  $cmd=$cmd."java -Xmx4g -jar $gatk_root/GenomeAnalysisTK.jar -T ApplyRecalibration  -R $genome -input $hc/$sample.vcf -mode SNP --ts_filter_level 99.0 -recalFile $vqsr/recalibrate_SNP.recal -tranchesFile $vqsr/recalibrate_SNP.tranches -o $vqsr/recalibrated_snps_raw_indels.vcf 2>&1 | tee $vqsr/ApplyRecalibration.SNP.log \n";
  $cmd=$cmd."java -Xmx4g -jar $gatk_root/GenomeAnalysisTK.jar -T VariantRecalibrator -R $genome -input $vqsr/recalibrated_snps_raw_indels.vcf -resource:mills,known=true,training=true,truth=true,prior=12.0 $mill_indel -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile $vqsr/recalibrate_INDEL.recal -tranchesFile $vqsr/recalibrate_INDEL.tranches -rscriptFile $vqsr/recalibrate_INDEL_plots.R\n";  
  $cmd=$cmd."java -Xmx4g -jar $gatk_root/GenomeAnalysisTK.jar -T ApplyRecalibration  -R $genome -input $vqsr/recalibrated_snps_raw_indels.vcf -mode INDEL --ts_filter_level 99.0 -recalFile $vqsr/recalibrate_INDEL.recal -tranchesFile $vqsr/recalibrate_INDEL.tranches -o $vqsr/$sample.filtered.vcf 2>&1 | tee $vqsr/ApplyRecalibration.INDEL.log\n";
  unless (-e $err) {create_fullshell ($shell, $qname, $err, $log, $cmd);  `qsub $shell >> ./jobs` ;}
  checkjobid("");  
}
sub vep {
  our $hc;
  our $perlbin_root;
  our $sample;;
  our $vep;  
  our $vep_root;
  our $vqsr;
  my $self = shift;	    
  my $cmd;
  my $err;
  my $qname = "vep";
  my $log;  
  my $set;
  my $shell;
  $cmd = "";
  $shell = "$vep/vep.sh";
  $err   = "$vep/vep.err";
  $log   = "$vep/vep.log";    	
  #$cmd=$cmd."$vep_root/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $vqsr/$sample.filtered.vcf --cache --everything --dir $vep_root/.vep/ --cache_version 81 --pick --vcf --fork 30 --offline --ASSEMBLY GRCh37 -o $vep/$sample-variant_effect_output.txt \n";
  $cmd=$cmd."$vep_root/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $hc/$sample.vcf --cache --everything --dir $vep_root/.vep/ --cache_version 81 --pick --vcf --fork 30 --offline --ASSEMBLY GRCh37 -o $vep/$sample-variant_effect_output.txt \n";
  $cmd=$cmd."$perlbin_root/vcf2txt.v2.pl $vep/$sample-variant_effect_output.txt > $vep/$sample-variant_effect_output.txt.1 \n";
  $cmd=$cmd."$perlbin_root/tab2xls.pl $vep/$sample-variant_effect_output.txt.1 $vep/$sample-variant_effect_output.xls \n";
  unless (-e $err) {create_fullshell ($shell, $qname, $err, $log, $cmd);  `qsub $shell >> ./jobs` ;} # comment on 2017/7/21 because there is an error when submit vep.sh
  #unless (-e $err) {create_fullshell ($shell, $qname, $err, $log, $cmd);  `sh $shell` ;}
  checkjobid("");
}
sub report {  
  our $bed;
  our $bedtools_root;
  our $bqsr;  
  our $dir; 
  our $fastqstats_root;
  our $fastqstats_template; 
  our $hc;
  our $index_template;
  our $perlbin_root;
  our $q30 ;
  our $report;  
  our $rscript;  
  our $sample;
  our $samstats_template;  
  our $stats_out;
  our $vep;   
  my $self = shift;   
  my $bam;  
  my $cmd;
  my $content;
  my $data;       # data structure for fastq-stats output     
  my $err;
  my $log;
  my $output;     # the dir of each set in report layer
  my $qname = "Q30";
  my @samples;   
  my $set;  
  my $shell;
  my $r1;
  my $r2;

  push (@samples, $sample); # 便宜行事, 這樣不用改 createFastqStatsHtml 和 createSamStatsHtml
  
# 原則上個別組的 data 放在各組的資料夾 report/sample
# 共有的 檔案放在 report/html
  $self -> inithtml ($report);    
  $output = "$report/$sample";	 # the dir of each set in report layer
  unless (-e "$output") {system ("mkdir -p $output");}  	
  # fastq statistic ; run fastq-stats of ea-utils to calc fq statistic use SGE
  $r1 = "$dir/reads/$sample\.R1.fastq.gz";
  $r2 = "$dir/reads/$sample\.R2.fastq.gz";	  	  	  
  $bam = "$bqsr/$sample.sorted.dedup.realign.recal.bam";
  fastqstats($fastqstats_root, $sample, $r1, $r2, $stats_out); 
  samstats($fastqstats_root, $sample, $bam, $stats_out);
  
  $shell = "$stats_out/$sample.q30.sh";
  $err   = "$stats_out/$sample.q30.err";
  $log   = "$stats_out/$sample.q30.log";
	
  $cmd=$cmd."python $q30 $r1 > $stats_out/$sample.R1.q30.log \n"; 
  $cmd=$cmd."python $q30 $r2 > $stats_out/$sample.R2.q30.log \n";
  unless (-e $err) {create_fullshell ($shell, $qname, $err, $log, $cmd);  `qsub $shell >> ./jobs` ;}
  checkjobid(""); 
      
  # process the output of fastq-stats	 & sam-stats!! Need to pass array reference to a subroutine of a module
  $data = procFastqStats (\@samples, $stats_out);			
  unless (-e "$output/fastqstats.html") {createFastqStatsHtml (\@samples, $data, $fastqstats_template, "$output");} # create html output of fastq-stats	
  $data = procSamStats (\@samples, $stats_out);			
  unless (-e "$output/samstats.html") {createSamStatsHtml (\@samples, $data, $samstats_template, "$output");} # create html output of sam-stats
	
  $content = $content."<li> $sample </a></li>";
  $content = $content."<li><a href=$sample/qc.html target=display>QC Summary</a></li>";
  $content = $content."<li><a href=$sample/fastqstats.html target=display>Fastq Stats</a></li>";
  $content = $content."<li><a href=$sample/samstats.html target=display>Alignment Stats</a></li>";
  
  unless (-e "$output/$sample-covreage.jpg") {targetcoverage ($sample, $bam, $rscript, $bedtools_root, $bed, $stats_out, $output, "wes");}
  $content = $content."<li><a href=$sample/$sample-covreage.jpg target=display>$sample coverage stat</a></li>";
  
  $self -> qstat;
  $self -> coverage;
  
  $self-> generateqc ($output);
  
  unless (-e "$output/$sample-variant_effect_output.xls") {system("cp $vep/$sample-variant_effect_output.xls $output/$sample-variant_effect_output.xls");} # copy variant_effect_output.xls file
  $content = $content."<li><a href=$sample/$sample-variant_effect_output.xls target=blank>Annotation</a></li>";

  unless (-e "$output/$sample-variant_effect_output.txt_summary.html") {system("cp $vep/$sample-variant_effect_output.txt_summary.html $output");} # copy variant_effect_output.txt_summary.html file
  $content = $content."<li><a href=$sample/$sample-variant_effect_output.txt_summary.html target=blank>Annotation Summary</a></li>";	    
  
  unless (-e "$output/$sample.vcf") {system("cp $hc/$sample.vcf $output");} # copy vcf
  $content = $content."<li><a href=$sample/$sample.vcf target=blank>VCF file</a></li>";
  
  $cmd = "sed 's|LeftList|$content|' $index_template > $report/index.html";
  process_cmd ($cmd);    
}
sub qstat {
  our $qualdata;
  our $sample;  
  our $stats_out;
  my $self = shift;  
  my @tags = ("R1", "R2");
  my $ln;
  my $log;
  my $rawbase = 0;
  my $tag;
  foreach my $tag (@tags) {
    $log = "$stats_out/$sample.$tag.q30.log";
	open (IN, '<', $log ) or die "Can't read $log file";
	while (<IN>) {
	  $ln = $_;
	  chomp($ln);
	  if ($ln=~/\('total bases:', (.*)\)/) {  
		$qualdata->{$sample}->{$tag}->{total_bases} = $1;
		$rawbase = $rawbase + $1;
	  }
	  if ($ln=~/\('q30 percents:', (.*)\)/) {
		$qualdata->{$sample}->{$tag}->{q30_percents} = $1;
		$rawbase = $rawbase + $1;
	  }		
	}
	close IN;
  }
  $qualdata->{$sample}->{raw_coverage} = $rawbase/45000000;
}
sub coverage {
  our $qualdata;
  our $sample;
  our $stats_out;
  my $self = shift;   
  my $above20;
  my $depth;  
  my $ln;
  my $hist;
  my $refbase;
  my $region;
  my $ratio;
  my $totalcoverbase;
  my @tab;
  
  $totalcoverbase = 0;
  $above20 = 0;
  $hist = "$stats_out/$sample.all.exome.coverage.hist.txt";
  open (IN, '<', $hist ) or die "Can't read $hist file";
  while (<IN>) {
	$ln = $_;
	chomp($ln);
	@tab = split(/\t/, $ln); 
	$depth = $tab[1];
	$region = $tab[2];
	$refbase = $tab[3];
	$ratio = $tab[-1];
	if ($depth >= 20) {
	  $above20 = $above20 +  $ratio;
	}
	$totalcoverbase = $totalcoverbase + ($depth * $region);
  }
  close IN;	
  $qualdata->{$sample}->{mean_coverage} = ($totalcoverbase/$refbase);
  $qualdata->{$sample}->{above20} = 100 * $above20;  
}
sub generateqc {
  our $qc_template;
  our $qualdata;
  our $sample;
  my $self = shift;
  my $output = shift;
  my $cmd;
  my $content;
  my $above20;
  my $mean_coverage;
  my @tags = ("R1", "R2");
  my $total_bases;
  my $q30_percents;  
  my $tag;
  my $raw_coverage;
  $content=$content."<table><tr><td>Tag</td><td>Total Base</td><td>% of Q30</td></tr>";
  foreach $tag (@tags) {
	$total_bases = $qualdata->{$sample}->{$tag}->{total_bases};
	$q30_percents = $qualdata->{$sample}->{$tag}->{q30_percents};
   $content=$content."<tr><td>$tag</td><td>".commify($total_bases)."</td><td>".commify($q30_percents)."</td></tr>";;
  }
  $content=$content."</table><br>";

  $content=$content."<table><tr><td>Raw Coverage</td><td>Mean coverage (X) - OnTarget</td><td>% of coverage at least 20x - OnTarget</td></tr>";
  $mean_coverage = $qualdata->{$sample}->{mean_coverage};
  $above20 = $qualdata->{$sample}->{above20};
  $raw_coverage = $qualdata->{$sample}->{raw_coverage};
  $content=$content."<tr><td>".commify($raw_coverage)."</td><td>".commify($mean_coverage)."</td><td>".commify($above20)."</td></tr>";;
  $content=$content."</table>";  
  
  
  $cmd = "sed 's|Content|$content|' $qc_template | sed 's|Title|QC|'> $output/qc.html";
  process_cmd ($cmd);
}
sub commify{
  my $number = shift;
  my $result;
  if ($number > 1000) {
     $number =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
  } else {
     $number = sprintf("%.2f", $number);
  }
  return $number;
}
sub inithtml {
  my $self = shift;
  my $output = shift;
  my $cmd;
  my @cmds;
  unless (-e "$output/images") {system("mkdir -p $output/images");} # image files
  unless (-e "$output/html")   {system("mkdir -p $output/html");} # html files	
  unless (-e "$output/files")  {system("mkdir $output/files");}   # files
  $cmd = "cp -r /export/EC1680U/Pipeline/Template/WES/css $output";
  push (@cmds, $cmd);
  $cmd = "cp -r /export/EC1680U/Pipeline/Template/WES/images/* $output/images";
  push (@cmds, $cmd);
  $cmd = "cp -r /export/EC1680U/Pipeline/Template/WES/DataAnalysisWorkflow.html $output/html/";
  push (@cmds, $cmd);
  $cmd = "cp -r /export/EC1680U/Pipeline/Template/WES/ExperimentWorkflow.html $output/html/";
  push (@cmds, $cmd);
  $cmd = "cp -r /export/EC1680U/Pipeline/Template/WES/help.pdf $output/files/";
  push (@cmds, $cmd);
  process_cmds_parallel(@cmds);
}
1;


