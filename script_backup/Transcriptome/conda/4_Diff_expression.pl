#! /usr/bin/env perl
# Author: kentchen
# Create date: 2017-06-12
# Last update: 2017-10-02

use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use Time::Seconds;


my $USER = `id -u -n`;
chomp($USER);
my ($input, $comp);
GetOptions(
    "input=s"  => \$input,
    "comp=s"    =>  \$comp,
);

if (! $comp){
    die "\n### ERROR: The comparison file is necessary\n";
} elsif (! $input){
    die "\n### ERROR: The input file is necessary\n";
} else {
    open OUT, ">>CMD_history.txt";
    print OUT CurrTime() . "$0 -i $input -c $comp\n";
    close OUT;
}


my $Trinity_home = "/usr/local/src/trinityrnaseq-Trinity-v2.4.0";
my $docker_basal_cmd = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` trinityrnaseq/trinityrnaseq:latest";
my $pwd = `pwd -L`;
chomp($pwd);
my $tab2xlsx = "/export/EC1680U/perl/bin/tab2xlsx.pl";


### Create outermost folder
my @folders = qw(1_AssemblyStats 2_Quantitation 3_DiffExpression);
my ($stat_dir, $quant_dir, $diff_dir) = @folders;
mkdir("$diff_dir");

### Generate metrics
print CurrTime() . "Generate metrix for genes and trans...\n";
my (@sample_trans, @sample_genes);
open my $sample_in, "<$input" or die "### ERROR: please check the '$input' existed or not\n";
while(<$sample_in>){
    chomp($_);
    my @sp = split/\t/, $_;
    my ($cond, $sample, $r1, $r2) = @sp;
    push(@sample_trans, "$quant_dir/$sample/RSEM.isoforms.results.plus1.tmp");
}
close $sample_in;

# trans matrix
my $expression_metrics_trans_cmd = "$docker_basal_cmd $Trinity_home/util/abundance_estimates_to_matrix.pl ".
"--est_method RSEM --name_sample_by_basedir --out_prefix $diff_dir/Trinity_trans @sample_trans";
system("$expression_metrics_trans_cmd");
system("sudo rm -rf $diff_dir/*TPM*");
system("$tab2xlsx $diff_dir/Trinity_trans.counts.matrix $diff_dir/Trinity_trans.counts.matrix.xlsx");
system("$tab2xlsx $diff_dir/Trinity_trans.TMM.EXPR.matrix $diff_dir/Trinity_trans.TMM.EXPR.matrix.xlsx");


### Differential expression main program
print CurrTime() . "Differential expression main program...\n";
open IN, "<$comp";
my ($comparisons, $cmd, @reports, @sig_reports);
my $raw_matrix = "$diff_dir/Trinity_trans.counts.matrix";
my $TMM_matrix = "$diff_dir/Trinity_trans.TMM.EXPR.matrix";
my $rscript = "/export/EC1680U/software/anaconda2/bin/Rscript";
while(<IN>){
    chomp;
    my @tab = split/\t/, $_;
    my $ctrl_count = scalar (split/,/, $tab[0]);
    my $case_count = scalar (split/,/, $tab[1]);
    my $ctrls = join '-', (split/,/, $tab[0]);
    my $cases = join '-', (split/,/, $tab[1]);
    $cmd = 'df <- read.table("'.$raw_matrix.'", row.names=1, header=T, sep="\t", check.names=F);';
#	$cmd .=  'df[df == 0] <- 0.01;'; # to avoid logFC that zero FPKM will be turned out to nothings
    $cmd .= 'ctrl.raw.df <- unlist(strsplit("'.$tab[0].'", "[,]"));';
    $cmd .= 'ctrl.raw.df <- df[, ctrl.raw.df, drop=FALSE];';
    $cmd .= 'mean.ctrl.raw.df <- as.data.frame(rowSums(ctrl.raw.df)/'.$ctrl_count.');';
    $cmd .= 'colnames(mean.ctrl.raw.df) <- c("'.$ctrls.'");';
    $cmd .= 'case.raw.df <- unlist(strsplit("'.$tab[1].'", "[,]"));';
    $cmd .= 'case.raw.df <- df[, case.raw.df, drop=FALSE];';
    $cmd .= 'mean.case.raw.df <- as.data.frame(rowSums(case.raw.df)/'.$case_count.');';
    $cmd .= 'colnames(mean.case.raw.df) <- c("'.$cases.'");';
    $cmd .= 'comp <- cbind(mean.ctrl.raw.df, mean.case.raw.df);';
#   $cmd .= 'redup <- SET2[, !duplicated((colnames(SET2)))];';
    $cmd .= 'write.table(comp, file="'."Trans.raw.counts.matrix".'", sep="\t", quote=F, row.names=T, col.names=NA);';
    system($rscript . " -e '" . $cmd . "'");

    $cmd = 'df <- read.table("'.$TMM_matrix.'", row.names=1, header=T, sep="\t", check.names=F);';
#	$cmd .=  'df[df == 0] <- 0.01;'; # to avoid logFC that zero FPKM will be turned out to nothings
    $cmd .= 'ctrl.tmm.df <- unlist(strsplit("'.$tab[0].'", "[,]"));';
    $cmd .= 'ctrl.tmm.df <- df[, ctrl.tmm.df, drop=FALSE];';
    $cmd .= 'mean.ctrl.tmm.df <- as.data.frame(rowSums(ctrl.tmm.df)/'.$ctrl_count.');';
    $cmd .= 'colnames(mean.ctrl.tmm.df) <- c("'.$ctrls.'");';
    $cmd .= 'case.tmm.df <- unlist(strsplit("'.$tab[1].'", "[,]"));';
    $cmd .= 'case.tmm.df <- df[, case.tmm.df, drop=FALSE];';
    $cmd .= 'mean.case.tmm.df <- as.data.frame(rowSums(case.tmm.df)/'.$case_count.');';
    $cmd .= 'colnames(mean.case.tmm.df) <- c("'.$cases.'");';
    $cmd .= 'comp <- cbind(mean.ctrl.tmm.df, mean.case.tmm.df);';
#   $cmd .= 'redup <- SET2[, !duplicated((colnames(SET2)))];';
    $cmd .= 'write.table(comp, file="'."Trans.TMM.counts.matrix".'", sep="\t", quote=F, row.names=T, col.names=NA);';
    system($rscript . " -e '" . $cmd . "'");
	
    my $comparisons = "$tab[0]-$tab[1]";
    mkdir("$diff_dir/$comparisons");

    system("mv Trans.*.counts.matrix $diff_dir/$comparisons");

	# run_DE_analysis.pl
	my $Trinity_home = "/usr/local/src/trinityrnaseq-Trinity-v2.4.0";
    my $docker_basal_cmd = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` trinityrnaseq/trinityrnaseq:latest";
    my $diff_edgeR_trans_cmd = "$docker_basal_cmd $Trinity_home/Analysis/DifferentialExpression/run_DE_analysis.pl ".
    "--matrix 3_DiffExpression/$comparisons/Trans.raw.counts.matrix --method edgeR ".
    "--output 3_DiffExpression/$comparisons --dispersion 0.1";
    system("$diff_edgeR_trans_cmd");

	# analyze_diff_expr.pl
    my $docker_mod_trans_cmd = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd`/$diff_dir/$comparisons trinityrnaseq/trinityrnaseq:latest";
    my $analyze_diff_trans_cmd = "$docker_mod_trans_cmd $Trinity_home/Analysis/DifferentialExpression/analyze_diff_expr.pl ".
    "--matrix Trans.TMM.counts.matrix ".
    "-P 0.05 -C 1";
    system("$analyze_diff_trans_cmd");
    system("sudo chown -R $USER:$USER $diff_dir");

	system("$tab2xlsx $diff_dir/$comparisons/Trans.raw.counts.matrix $diff_dir/$comparisons/Trans.raw.counts.matrix.xlsx");
	system("$tab2xlsx $diff_dir/$comparisons/Trans.TMM.counts.matrix $diff_dir/$comparisons/Trans.TMM.counts.matrix.xlsx");

	# Generate up/down regulated transcripts list for following KEGG pathway
	system("cat $diff_dir/$comparisons/*.$ctrls-UP.subset | cut -f1 | tail -n +2 > $diff_dir/$comparisons/$comparisons\_DOWN.list");
	system("cat $diff_dir/$comparisons/*.$cases-UP.subset | cut -f1 | tail -n +2 > $diff_dir/$comparisons/$comparisons\_UP.list");


	# Generate each comparison DE report
	my $de_result = `ls $diff_dir/$comparisons/*.DE_results`;
    chomp($de_result);
	open my $DE_file, "<$de_result";
	open my $DE_report, ">$diff_dir/$comparisons/$comparisons.DE_result";
	while(my $DE_txt = <$DE_file>){
		chomp($DE_txt);
		$DE_txt =~ s/sampleA/transcript_id/g;
		$DE_txt =~ s/logCPM/logFC ($comparisons)/g;
		$DE_txt =~ s/FDR/pvalue ($comparisons)/g;
		my @DE_sp = split/\t/, $DE_txt;
		print $DE_report "$DE_sp[0]\t$DE_sp[3]\t$DE_sp[5]\n";
	}
	push(@reports, "$diff_dir/$comparisons/$comparisons.DE_result");
	close $DE_file;
	close $DE_report;
	
	system("$tab2xlsx $de_result $de_result.xlsx");

	# Generate each comparison P<0.05 & logFC>1 or logFC<-1 DE report for grouped heatmap plot
	my $de_sig_result = `ls $diff_dir/$comparisons/*DE_results.P0.05_C1.DE.subset`;
	chomp($de_sig_result);
	open my $DE_sig_file, "<$de_sig_result";
	open my $DE_sig_report, ">$diff_dir/$comparisons/$comparisons.P0.05_C1.DE_result";
	while(my $DE_sig_txt = <$DE_sig_file>){
		chomp($DE_sig_txt);
        $DE_sig_txt =~ s/sampleA/transcript_id/g;
        $DE_sig_txt =~ s/logCPM/logFC ($comparisons)/g;
		my @DE_sig_sp = split/\t/, $DE_sig_txt;
		print $DE_sig_report "$DE_sig_sp[0]\t$DE_sig_sp[3]\n";
	}
	push(@sig_reports, "$diff_dir/$comparisons/$comparisons.P0.05_C1.DE_result");
	close $DE_sig_file;
	close $DE_sig_report;
	
	my @subsets = `ls $diff_dir/$comparisons/*.subset`;
	chomp(@subsets);
	foreach my $subset (@subsets){
		system("$tab2xlsx $subset $subset.xlsx");
	}

	# Generate Scatter Plot (all & significant only)
	system("cut -f2,3 $diff_dir/$comparisons/Trans.TMM.counts.matrix > $diff_dir/$comparisons/Trans.TMM.counts.mod.matrix");
	my $scatter_de_result = `ls $diff_dir/$comparisons/*.edgeR.DE_results`; chomp($scatter_de_result);
	my $scatter_sig_de_result = `ls $diff_dir/$comparisons/*.DE.subset`; chomp($scatter_sig_de_result);
	system("paste -d\"\t\" $scatter_de_result $diff_dir/$comparisons/Trans.TMM.counts.mod.matrix > $diff_dir/$comparisons/Trans.DE.forScatter.txt");
	system("cp /export/EC1680U/Rscript/Trinity/Scatter.Trinity.r .");
	$tab[0] =~ s/,/-/g;
	system("sed -i 's/sample1/$tab[0]/g' Scatter.Trinity.r");
	$tab[1] =~ s/,/-/g;
	system("sed -i 's/sample2/$tab[1]/g' Scatter.Trinity.r");
	system("$rscript Scatter.Trinity.r $diff_dir/$comparisons/Trans.DE.forScatter.txt scatter_all.html");
	system("$rscript Scatter.Trinity.r $scatter_sig_de_result scatter_sig.html");
	system("rm -rf $diff_dir/$comparisons/scatter*");
	system("mv -f scatter_all* scatter_sig* $diff_dir/$comparisons"); # savewidget could not output to specific folder
}
close IN;

### Generate group comparison table
print CurrTime() . "Generate group comparison table...\n";
my $merge_logFC_cmd;
my @logFC_list;
my $num = 1;
foreach my $sig_report (@sig_reports){
	my $name = "tmp$num";
	$merge_logFC_cmd .= $name.' <- read.table("'.$sig_report.'", header=T, sep="\t", check.names=FALSE, fill=TRUE, quote="\"");';
	push(@logFC_list, $name);
	$num++;
}
my $logFC_list_comma = join(',', @logFC_list);
$merge_logFC_cmd .= 'merge <- Reduce(function(x,y) merge(x, y, by="transcript_id", all=F), list('.$logFC_list_comma.'));';
$merge_logFC_cmd .= 'write.table(merge, file="'.$diff_dir.'/allGroups_comparison.report", sep="\t", quote=F, row.names = FALSE);';
system($rscript . " -e '" . $merge_logFC_cmd . "'");
system("sed -i 's/logFC (//g' $diff_dir/allGroups_comparison.report");
system("sed -i 's/)//g' $diff_dir/allGroups_comparison.report");
system("$tab2xlsx $diff_dir/allGroups_comparison.report $diff_dir/allGroups_comparison.report.xlsx");

### Plot inter-sample heatmap
print CurrTime() . "Plot inter-sample heatmap intersection...\n";
if(scalar @sig_reports >= 2){
	system("$rscript /export/EC1680U/Rscript/Trinity/heatmap_eachSample.R $diff_dir/Trinity_trans.TMM.EXPR.matrix $diff_dir");
}
else {
	system("touch $diff_dir/noSample_heatmap_plot");
}

### Plot PCA
print CurrTime() . "Plot PCA...\n";
my $PCA_cmd = "/usr/bin/Rscript /export/EC1680U/Rscript/Trinity/PCA_plot.R ".
"$diff_dir/Trinity_trans.TMM.EXPR.matrix $diff_dir/all_samples_PCA.png";
system("$PCA_cmd");

# Plot inter-group heatmap
print CurrTime() . "Plot inter-group heatmap...\n";
if(scalar @sig_reports >= 2){
	system("$rscript /export/EC1680U/Rscript/Trinity/heatmap_eachGroup.R $diff_dir/allGroups_comparison.report $diff_dir");
}
else {
	system("touch $diff_dir/noGroup_heatmap_plot");
}

### Plot the correlation within all samples
print CurrTime() . "Plot correlation within all samples...\n";
if(scalar @sig_reports >= 2){
	system("$rscript /export/EC1680U/Rscript/Trinity/correlation_allSample.R $diff_dir/Trinity_trans.TMM.EXPR.matrix $diff_dir/all_samples_correlation.png");
}
else {
    system("touch $diff_dir/noSample_correlation_plot");
}



### Final DE report
print CurrTime() . "Generate final DE report...\n";
my $merge_cmd;
my @list;
$num = 1;
foreach my $report (@reports){
	my $name = "tmp$num";
	$merge_cmd .= $name.' <- read.table("'.$report.'", header=T, sep="\t", check.names=FALSE, fill=TRUE, quote="\"");';
	push(@list, $name);
	$num++;
}
my $list_comma = join(',', @list);
$merge_cmd .= 'merge <- Reduce(function(x,y) merge(x, y, by="transcript_id", all=TRUE), list('.$list_comma.'));';
$merge_cmd .= 'merge <- merge[, !duplicated(colnames(merge))];';
$merge_cmd .= 'merge$length.y <- NULL;';
$merge_cmd .= 'write.table(merge, file="'.$diff_dir.'/all_DE_comparison.report.tmp", sep="\t", quote=F, na="-", row.names = FALSE);';
system($rscript . " -e '" . $merge_cmd . "'");

# Grouped logFC & PValue respectively
my $num_report = scalar @reports;
open my $all_DE_report, "<$diff_dir/all_DE_comparison.report.tmp";
open my $final_DE_report, ">$diff_dir/final_DE_comparison.report.txt";
while(<$all_DE_report>){
	chomp;
	my @sp = split/\t/, $_;
	my $interval = ($#sp+1)/$num_report;
	print $final_DE_report "$sp[0]\t";
	for(my $i=1; $i<=$#sp; $i+=$interval){
		print $final_DE_report "$sp[$i]\t";
	}
	for(my $j=2; $j<=$#sp+1; $j+=$interval){
		print $final_DE_report "$sp[$j]\t";
	}
	print $final_DE_report "\n";
}
close $all_DE_report;
close $final_DE_report;

# Remove un-use files in comparisons
open IN, "<$comp";
while(<IN>){
    chomp;
    my @tab = split/\t/, $_;
    my $ctrl_count = scalar (split/,/, $tab[0]);
    my $case_count = scalar (split/,/, $tab[1]);
    my $ctrls = join '-', (split/,/, $tab[0]);
    my $cases = join '-', (split/,/, $tab[1]);

	my $comparisons = "$tab[0]-$tab[1]";

	system("mkdir $diff_dir/$comparisons/tmp");
    system("mv $diff_dir/$comparisons/*.pdf $diff_dir/$comparisons/*.xlsx $diff_dir/$comparisons/Trans.*.counts.matrix $diff_dir/$comparisons/*.list ".
    "$diff_dir/$comparisons/*.DE_results $diff_dir/$comparisons/*.subset $diff_dir/$comparisons/*edgeR.count_matrix $diff_dir/$comparisons/*.html ".
    "$diff_dir/$comparisons/tmp/");
    system("ls $diff_dir/$comparisons/* | grep -v tmp | xargs rm -f");
    system("mv $diff_dir/$comparisons/tmp/* $diff_dir/$comparisons ; rm -rf $diff_dir/$comparisons/tmp");	
}
close IN;

system("sed -i 's/\t\$//g' $diff_dir/final_DE_comparison.report.txt");
system("$tab2xlsx $diff_dir/final_DE_comparison.report.txt $diff_dir/final_DE_comparison.report.xlsx");

system("rm -rf $diff_dir/all_DE_comparison.report.tmp");
system("sudo rm -rf __tmp_runTMM.R");



sub CurrTime {
    my $t = localtime;
    my $time = "[".$t->hms."]*** ";
    return "$time";
}
