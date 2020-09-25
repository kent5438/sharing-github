#! /usr/bin/env perl

use warnings;
use strict;

### Generate comparison for loop in rmarkdown
my $comp = $ARGV[0];
my @replaces;

open my $rmd, "<quantitation_template.Rmd";
open my $rmd_out, ">quantitation.Rmd";
while(my $line = <$rmd>){
	chomp($line);
	if($line =~ /###COMPARISON1###/){
		open IN, "<$comp";
		while(<IN>){
			chomp;
			my @tab = split/\t/, $_;
			my $ctrl_count = scalar (split/,/, $tab[0]);
			my $case_count = scalar (split/,/, $tab[1]);
			my $ctrls = join '-', (split/,/, $tab[0]);
			my $cases = join '-', (split/,/, $tab[1]);
			my $comparisons = "$tab[0]-$tab[1]";
			
            my $replace = "#### $comparisons\n\n";
            $replace .= "[All]\n\n";
			$replace .= "---> Beware! over 100k transcripts may cause PC crashing of this link! <---\n\n";
			$replace .= "<a href=\"files/3_DiffExpression/$comparisons/scatter_all.html\" target=\"_blank\">$comparisons: All (.html)</a>\n\n";
#            $replace .= '```{r, out.width = "700px", out.height="800px", echo=FALSE, message=FALSE, warning=FALSE}'."\n";
#            $replace .= "knitr::include_url(\"files/3_DiffExpression/$comparisons/scatter_all.html\")\n";
#            $replace .= '```'."\n";
            $replace .= "\n<br>\n\n";
            $replace .= "[Significant (FDR<0.05 / log2FC+-1)]\n\n";
			$replace .= "<a href=\"files/3_DiffExpression/$comparisons/scatter_sig.html\" target=\"_blank\">$comparisons: Sig (.html)</a>\n\n";
#            $replace .= '```{r, out.width = "700px", out.height="800px", echo=FALSE, message=FALSE, warning=FALSE}'."\n";
#            $replace .= "knitr::include_url(\"files/3_DiffExpression/$comparisons/scatter_sig.html\")\n";
#            $replace .= '```'."\n";
            $replace .= "\n<br>\n\n";
			push(@replaces,$replace);
		}
		close IN;
	}
	
	$line =~ s/###COMPARISON1###/@replaces/g;
	
	@replaces='';

	if($line =~ /###COMPARISON2###/){
		open IN, "<$comp";
		while(<IN>){
			chomp;
			my @tab = split/\t/, $_;
			my $ctrl_count = scalar (split/,/, $tab[0]);
			my $case_count = scalar (split/,/, $tab[1]);
			my $ctrls = join '-', (split/,/, $tab[0]);
			my $cases = join '-', (split/,/, $tab[1]);
			my $comparisons = "$tab[0]-$tab[1]";

			my $replace = "#### $comparisons\n\n";
			$replace .= "[MA-Plot (FDR<0.05)]\n\n";
			$replace .= "<a href=\"files/3_DiffExpression/$comparisons/glimma_MA/MA-plot.html\" target=\"_blank\">$comparisons: MA-Plot (.html)</a>\n\n";
			$replace .= "\n<br>\n\n";
			$replace .= "[Volcano-Plot (FDR<0.05)]\n\n";
			$replace .= "<a href=\"files/3_DiffExpression/$comparisons/glimma_Volcano/Volcano-plot.html\" target=\"_blank\">$comparisons: Volcano-Plot (.html)</a>\n\n";
			$replace .= "\n<br>\n\n";
			push(@replaces,$replace);
		}
		close IN;
	}
	$line =~ s/###COMPARISON2###/@replaces/g;

	print $rmd_out "$line\n";
}
close $rmd;
close $rmd_out;

system("sed -i 's/^ //g' quantitation.Rmd");
