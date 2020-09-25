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
	if($line =~ /###COMPARISON###/){
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
			$replace .= "<a href=\"files/3_DiffExpression/$comparisons/scatter_all.html\" target=\"_blank\">$comparisons: All (.html)</a>\n\n";
#            $replace .= '```{r, out.width = "700px", out.height="800px", echo=FALSE, message=FALSE, warning=FALSE}'."\n";
#            $replace .= "knitr::include_url(\"files/3_DiffExpression/$comparisons/scatter_all.html\")\n";
#            $replace .= '```'."\n";
            $replace .= "\n<br>\n\n";
            $replace .= "[Significant (pvalue<0.05 / FC +- 1)]\n\n";
			$replace .= "<a href=\"files/3_DiffExpression/$comparisons/scatter_sig.html\" target=\"_blank\">$comparisons: Sig (.html)</a>\n\n";
#            $replace .= '```{r, out.width = "700px", out.height="800px", echo=FALSE, message=FALSE, warning=FALSE}'."\n";
#            $replace .= "knitr::include_url(\"files/3_DiffExpression/$comparisons/scatter_sig.html\")\n";
#            $replace .= '```'."\n";
            $replace .= "\n<br>\n\n";
			push(@replaces,$replace);
		}
		close IN;
	}
	$line =~ s/###COMPARISON###/@replaces/g;
	print $rmd_out "$line\n";
}
close $rmd;
close $rmd_out;

system("sed -i 's/^ //g' quantitation.Rmd");
