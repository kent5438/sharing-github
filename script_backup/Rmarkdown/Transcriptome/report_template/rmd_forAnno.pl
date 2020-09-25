#! /usr/bin/env perl

use warnings;
use strict;

### Modified Transdecoder README.txt
system("rsync -avq /export/EC1680U/Pipeline/Trinity/report_template/README.txt files/4_Annotation/Gene_prediction");
my $count1 = `grep -c '>' files/1_AssemblyStats/Trinity.fasta`; chomp($count1);
my $count2 = `grep -c '>' files/1_AssemblyStats/Trinity.95.fasta`; chomp($count2);
my $count3 = `grep -c '>' files/4_Annotation/Gene_prediction/Trinity.95.fasta.transdecoder.cds`; chomp($count3);
system("sed -i 's/#COUNT1#/$count1/' files/4_Annotation/Gene_prediction/README.txt");
system("sed -i 's/#COUNT2#/$count2/' files/4_Annotation/Gene_prediction/README.txt");
system("sed -i 's/#COUNT3#/$count3/' files/4_Annotation/Gene_prediction/README.txt");

### link Trinity.fasta & Trinity.95.fasta to Gene_prediction
chdir("files/4_Annotation/Gene_prediction");
if((! -e "Trinity.fasta") && (! -e "Trinity.95.fasta")){
	system("ln -s ../../1_AssemblyStats/Trinity.fasta .");
	system("ln -s ../../1_AssemblyStats/Trinity.95.fasta .");
}
chdir("../../../");

### ### Generate comparison for loop in rmarkdown
my $comp = $ARGV[0];

my @replaces_GO;
my @replaces_KEGG;
open my $rmd, "<annotation_template.Rmd";
open my $rmd_out, ">annotation.Rmd";
while(my $line = <$rmd>){
    chomp($line);
    if($line =~ /###GO_COMPARISON###/){
        open IN, "<$comp";
        while(<IN>){
            chomp;
            my @tab = split/\t/, $_;
            my $ctrl_count = scalar (split/,/, $tab[0]);
            my $case_count = scalar (split/,/, $tab[1]);
            my $ctrls = join '-', (split/,/, $tab[0]);
            my $cases = join '-', (split/,/, $tab[1]);
            my $comparisons = "$tab[0]-$tab[1]";

			my $go_enrich_subset_all = `ls files/4_Annotation/GO/$comparisons\_enrich/*\.DE\.subset\.GOseq\.enriched\.xlsx`;
			chomp($go_enrich_subset_all);
			my $go_enrich_plot_all = `ls files/4_Annotation/GO/$comparisons\_enrich/*\.DE.subset\.GOseq\.enriched\.GOplot_dat\.pdf`;
			chomp($go_enrich_plot_all);
			my $go_enrich_subset_up = `ls files/4_Annotation/GO/$comparisons\_enrich/*\.$tab[1]-UP\.subset\.GOseq\.enriched\.xlsx`;
			chomp($go_enrich_subset_up);
			my $go_enrich_subset_down = `ls files/4_Annotation/GO/$comparisons\_enrich/*\.$tab[0]-UP\.subset\.GOseq\.enriched\.xlsx`;
            chomp($go_enrich_subset_down);

            my $replace = "#### $comparisons\n\n";
            $replace .= "[All regulation]\n\n";
            $replace .= "[$comparisons (.xlsx)]($go_enrich_subset_all)\n\n";
			$replace .= "<a href=\"$go_enrich_plot_all\" target=\"_blank\">GO enrichment bar plot (.pdf)</a>\n\n";
            $replace .= "***\n\n";
            
			$replace .= "[Sig UP (pvalue<0.05 / logFC >= 1)]\n\n";
            $replace .= "[$comparisons: UP (.xlsx)]($go_enrich_subset_up)\n\n";
            $replace .= "<br>\n\n";

			$replace .= "[Sig DOWN (pvalue<0.05 / logFC <= -1)]\n\n";
            $replace .= "[$comparisons: DOWN (.xlsx)]($go_enrich_subset_down)\n\n";
            $replace .= "<br>\n\n";

            push(@replaces_GO,$replace);
        }
        close IN;
    }
	$line =~ s/###GO_COMPARISON###/@replaces_GO/g;

    if($line =~ /###KEGG_COMPARISON###/){
        open IN, "<$comp";
        while(<IN>){
            chomp;
            my @tab = split/\t/, $_;
            my $ctrl_count = scalar (split/,/, $tab[0]);
            my $case_count = scalar (split/,/, $tab[1]);
            my $ctrls = join '-', (split/,/, $tab[0]);
            my $cases = join '-', (split/,/, $tab[1]);
            my $comparisons = "$tab[0]-$tab[1]";

            my $kegg_subset_up = `ls files/4_Annotation/KEGG/$comparisons\_kegg/ec2kegg_UP.xls`;
            chomp($kegg_subset_up);
            my $kegg_subset_down = `ls files/4_Annotation/KEGG/$comparisons\_kegg/ec2kegg_DOWN.xls`;
            chomp($kegg_subset_down);
            my $kegg_subset_up_plot = `ls files/4_Annotation/KEGG/$comparisons\_kegg/kegg_barchart_UP.png`;
            chomp($kegg_subset_up_plot);
			my $kegg_subset_down_plot = `ls files/4_Annotation/KEGG/$comparisons\_kegg/kegg_barchart_DOWN.png`;
			chomp($kegg_subset_down_plot);

            my $replace = "#### $comparisons\n\n";
            $replace .= "[Sig UP (pvalue<0.05 / logFC >= 1)]\n\n";
            $replace .= "[$comparisons: UP (.xls)]($kegg_subset_up)\n\n";
			$replace .= "```{r, out.width = \"600px\", echo=FALSE, message=FALSE, warning=FALSE}\n";
			$replace .= "knitr::include_graphics(\"$kegg_subset_up_plot\")\n";
			$replace .= "```\n\n";
            $replace .= "<br>\n\n";

            $replace .= "[Sig DOWN (pvalue<0.05 / logFC <= -1)]\n\n";
            $replace .= "[$comparisons: DOWN (.xlsx)]($kegg_subset_down)\n\n";
			$replace .= "```{r, out.width = \"600px\", echo=FALSE, message=FALSE, warning=FALSE}\n";
			$replace .= "knitr::include_graphics(\"$kegg_subset_down_plot\")\n";
			$replace .= "```\n\n";
            $replace .= "<br>\n\n";

            push(@replaces_KEGG,$replace);
        }
        close IN;
    }
    $line =~ s/###KEGG_COMPARISON###/@replaces_KEGG/g;

    print $rmd_out "$line\n";
}
close $rmd;
close $rmd_out;

system("sed -i 's/^ //g' annotation.Rmd");
