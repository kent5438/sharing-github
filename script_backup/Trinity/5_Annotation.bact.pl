#! /usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use Time::Seconds;
use lib '/export/EC1680U/kentchen/lib';
use Time qw(CurrTime);


my $pwd = `pwd -L`; chomp($pwd);
my $USER = `id -u -n`;
chomp($USER);
my $tab2xls = "/export/EC1680U/perl/bin/tab2xls.pl";
my $tab2xlsx = "/export/EC1680U/perl/bin/tab2xlsx.pl";

my ($help, $prefix, $code);
GetOptions(
	"help!"	=> \$help,
	"prefix=s"	=> \$prefix,
	"kegg=s"	=> \$code,
);

if($help){
	die "\n### USAGE ###
* Used for Pacbio assembled genome annotation pipeline by following tools
- Prokka (prokka)
- COG (rpsblast & cdd2cog.pl)
- GO (interproscan)
- KEGG (ec2kegg.pl)

E.g.)
% $0 -p PB17042-Y2U -k bsu

--prefix/-p : name of prefix
--kegg/-k : KEGG_organism_code (http://www.genome.jp/kegg/catalog/org_list.html)
"
} elsif (! $prefix){
	die "\n### ERROR: The output prefix is necessary, e.g. PB17042-Y2U\n";
} elsif (! $code){
	die "\n### ERROR: KEGG organism code is necessary, e.g. bsu\n";
} else {
	open OUT, ">>CMD_history.txt";
	print OUT CurrTime() . "$0 -p $prefix -k $code\n";
	close OUT;
}

### Create folders ###
print CurrTime() . "Create folder in $pwd...\n";
my @folders = qw(1_Assembly 2_Subreads 3_Prokka 4_Annotation 5_Report);
my ($assembly_dir, $subread_dir, $prokka_dir, $anno_dir, $report_dir) = @folders;
foreach my $folder (@folders){
	if (! -e "$folder"){mkdir "$folder";}
}
my $go_dir = "$anno_dir/GO";
my $kegg_dir = "$anno_dir/KEGG";
my $cog_dir = "$anno_dir/COG";



### Retrieve eC number from gff, used for the following ec2kegg ###
my $faa = "$prokka_dir/$prefix.faa";
my $ffn = "$prokka_dir/$prefix.ffn";
my $fna = "$prokka_dir/$prefix.fna";

print CurrTime() . "Retrieve eC number from gff...\n";
my $ec_cmd = "grep 'eC_number=' $prokka_dir/$prefix.gff | cut -f9 | awk -F';' '{print \$3}' | sed 's/eC_number=//g' > $prokka_dir/$prefix.ec";
system("$ec_cmd");

my $gene_cmd = "grep 'eC_number=' $prokka_dir/$prefix.gff | cut -f9 | awk -F';' '{print \$1}' | sed 's/ID=//g' > $prokka_dir/$prefix.geneID.tmp";
my $paste_cmd = "paste $prokka_dir/$prefix.geneID.tmp $prokka_dir/$prefix.ec > $go_dir/GeneID_ec.txt";
mkdir("$go_dir");
system("$gene_cmd; $paste_cmd");

### GO annotation work by Interproscan ###
print CurrTime() . "Start GO annotation work...\n";

# Interproscan work (we just want to parse GO data)
print CurrTime() . "Start Interproscan work... (SLOW!)\n";
my $interpro_docker_basal = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` kentchendocker/interproscan:latest";
my $interpro_cmd = "$interpro_docker_basal /interproscan-5.23-62.0/interproscan.sh ".
"-i $prokka_dir/$prefix.faa -appl Pfam -goterms -pa -b $prefix\_interpro -dp";	
if(! -e "$go_dir/GeneID_GO.txt"){
	system("$interpro_cmd");
	system("sudo chown -R $USER:users temp $prefix\_interpro.*");
	system("awk -F\"\t\" '{print \$1 \"\t\" \$14}' $prefix\_interpro.tsv > GeneID_GO.txt.tmp");
	my $cmd = 'df <- read.table("GeneID_GO.txt.tmp", header=F, sep="\t");';
	$cmd .= 'df <- df[!(is.na(df$V2) | df$V2==""), ];';
	$cmd .= 'a <- aggregate(V2~V1, data=df, paste, collapse="|");';
	$cmd .= 'write.table(a, "'.$go_dir.'/GeneID_GO.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE);';
	system("Rscript" . " -e '" . $cmd . "'");
	system("sed -i 's/|/,/g' $go_dir/GeneID_GO.txt");
#	system("sed -i -re 's/_[0-9]+\t/\t/g' $go_dir/GeneID_GO.txt");
	system("sed -i '1i Prokka_ID\tGOs' $go_dir/GeneID_GO.txt");
	system("rm -rf *_interpro.* GeneID_GO.txt.tmp temp");
}


# Parse GO slim (generic) by Trinotate
my $go_slim_mapping = "Trinotate_GO_to_SLIM.pl $go_dir/GeneID_GO.txt > $go_dir/GO_slim_mapping.txt";
if(! -e "$go_dir/GO_slim_mapping.txt"){system("$go_slim_mapping");}

open IN, "<$go_dir/GO_slim_mapping.txt";
open OUT, ">$go_dir/GO_mapping.txt.tmp";
while(<IN>){
    chomp;
    my @sp = split/\t/, $_;
    no warnings 'uninitialized';
    if(($sp[2] ne "biological_process") && ($sp[2] ne "cellular_component") && ($sp[2] ne "molecular_function")){
        print OUT "$sp[0]\t$sp[1]\t$sp[2]\t$sp[3]\t$sp[4]\n";
    }
}
close IN;
close OUT;

# remove 0 hit of GO_mapping.txt
my $cmd = 'df <- read.table("'.$go_dir.'/GO_mapping.txt.tmp", header=F, sep="\t", fill=TRUE, quote="\"");';
$cmd .= 'df <- df[!(is.na(df$V4) | df$V4=="0"), ];';
$cmd .= 'write.table(df, "'.$go_dir.'/GO_mapping.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE);';
system("Rscript" . " -e '" . $cmd . "'");

system("sed -i '1i Class\tGO_ID\tFunction\tCount\tDescription' $go_dir/GO_mapping.txt");

# Plot GO barchart
print CurrTime() . "Plot GO mapping graph...\n";
my $GO_mapping_R = "Rscript /export/EC1680U/Rscript/Pacbio/go_mapping_bar_chart.R $go_dir/GO_mapping.txt $go_dir &> /dev/null";
if(! -e "$go_dir/GO_barcahrt.png"){system("$GO_mapping_R");}
system("rm -rf $go_dir/GO_slim_mapping.txt $go_dir/GO_mapping.txt.tmp");

my @go_txt = `ls $go_dir/*.txt`;
chomp(@go_txt);
foreach my $each_go_txt (@go_txt){
    system("$tab2xlsx $each_go_txt $each_go_txt.xlsx");
}


### ec2kegg annotation work ###
if(! -e $kegg_dir){mkdir("$kegg_dir");}
my $ec2kegg_bin = "/export/EC1680U/software/ec2kegg";
my $ec2kegg_update_cmd = "$ec2kegg_bin/get_definitions.pl";
my $ec2kegg_cmd = "$ec2kegg_bin/ec2kegg.pl $code $prokka_dir/$prefix.ec > $kegg_dir/$prefix\_ec2kegg.txt";

print CurrTime() . "Update KEGG definition database...\n";
system("$ec2kegg_update_cmd");

print CurrTime() . "Parse KEGG annotation from eC number...\n";
system("$ec2kegg_cmd");
system("sed 's/||/\t/g' $kegg_dir/$prefix\_ec2kegg.txt > $kegg_dir/$prefix\_ec2kegg.txt.tmp");
system("grep -v '^00000' $kegg_dir/$prefix\_ec2kegg.txt.tmp > $kegg_dir/$prefix\_ec2kegg.txt");
system("rm -rf $kegg_dir/$prefix\_ec2kegg.txt.tmp");
system("$tab2xls $kegg_dir/$prefix\_ec2kegg.txt $kegg_dir/$prefix\_ec2kegg.xls");

print CurrTime() . "Plot KEGG bar chart...\n";
my $plot_kegg_cmd = "Rscript /export/EC1680U/Rscript/Pacbio/kegg_horizontal_bar_chart.R ".
"$kegg_dir/$prefix\_ec2kegg.txt $kegg_dir $code &> /dev/null";
system("$plot_kegg_cmd");

system("cp /export/EC1680U/kentchen/Doc/ec2kegg_README.txt $kegg_dir");

# Add pathway name into GeneID_pathway.txt
open my $pathway_out, ">$kegg_dir/GeneID_pathway.txt.tmp";
my $ec_file = "$go_dir/GeneID_ec.txt";
my @enzymes = `cut -f2 $ec_file | sed 's/\.-//g'`;
chomp(@enzymes);

system("cut -f2,3,13 $kegg_dir/$prefix\_ec2kegg.txt > $kegg_dir/$prefix\_ec2kegg.txt.tmp");

my $ec2kegg_file = "$kegg_dir/$prefix\_ec2kegg.txt.tmp";
foreach my $enzyme (@enzymes){
	my $output = `egrep '$enzyme,|,$enzyme|,$enzyme,' $ec2kegg_file | cut -f1 | tr "\n" , | sed 's/,\$//'`;
	chomp($output);
	print $pathway_out "$output\n";
}
close $pathway_out;

system("paste $go_dir/GeneID_ec.txt $kegg_dir/GeneID_pathway.txt.tmp | cut -f1,3 > $kegg_dir/GeneID_pathway.txt");
system("sed -i '1i Prokka_ID\tPathway' $kegg_dir/GeneID_pathway.txt");
system("rm -f $kegg_dir/*.txt.tmp");


### COG annotation work ###
if(! -e $cog_dir){mkdir("$cog_dir");}
print CurrTime() . "COG annotation work... [rpsbalst stage]\n";
my $cog_db = "/export/EC1680U/DataBase/COG";
my $rpsblast_cmd = "rpsblast -query $prokka_dir/$prefix.faa ".
"-db $cog_db/Cog -out $cog_dir/$prefix\_rpsblast.out ".
"-evalue 1e-2 -outfmt 6 -num_threads 32";
if(! -e "$cog_dir/$prefix\_rpsblast.out"){system("$rpsblast_cmd");}

print CurrTime() . "COG annotation work... [Parse COG stage]\n";
my $cdd2cog = "$cog_db/cdd2cog.pl";
my $cdd2cog_cmd = "perl $cdd2cog -r $cog_dir/$prefix\_rpsblast.out ".
"-c $cog_db/cddid.tbl -f $cog_db/fun.txt -w $cog_db/whog -a";
system("$cdd2cog_cmd");
system("mv results/* $cog_dir");
system("rm -rf results");

open IN, "<$cog_dir/func_stats.txt";
open OUT, ">$cog_dir/func_stats.txt.tmp";
while(<IN>){
    chomp;
    my @sp = split/\t/, $_;
    print OUT "$sp[0]\t$sp[0]: $sp[1]\t$sp[2]\n";
}
close IN;
close OUT;
system("rm -rf $cog_dir/func_stats.txt");
system("mv $cog_dir/func_stats.txt.tmp $cog_dir/func_stats.txt");

print CurrTime() . "COG annotation work... [Plot COG bar chart...]\n";
my $plot_cog_cmd = "Rscript /export/EC1680U/Rscript/Pacbio/cog_horizontal_bar_chart.R ".
"$cog_dir/func_stats.txt $cog_dir &> /dev/null";
if(! -e "$cog_dir/cog_barchart.png"){
    system("$plot_cog_cmd");
}

if (! -e "$cog_dir/GeneID_COG.txt"){
	system("awk -F\"\t\" '{print \$1 \"\t\" \$2}' $cog_dir/protein-id_cog.txt > $cog_dir/protein-id_cog.txt.tmp");
	my $cmd = 'df <- read.table("'.$cog_dir.'/protein-id_cog.txt.tmp", header=F, sep="\t");';
	$cmd .= 'df <- df[!(is.na(df$V2) | df$V2==""), ];';
	$cmd .= 'a <- aggregate(V2~V1, data=df, paste, collapse=",");';
	$cmd .= 'write.table(a, "'.$cog_dir.'/GeneID_COG.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE);';
	system("Rscript" . " -e '" . $cmd . "'");
	system("sed -i '1i Prokka_ID\tCOGs' $cog_dir/GeneID_COG.txt");
	system("rm -rf $cog_dir/protein-id_cog.txt.tmp");
}

my @cog_txt = `ls $cog_dir/*.txt`;
chomp(@cog_txt);
foreach my $each_cog_txt (@cog_txt){
    system("$tab2xlsx $each_cog_txt $each_cog_txt.xlsx");
}

## Using BLASTP to fill up the predicted hypothetical protein from prokka
print CurrTime() . "Using BLASTP to fill up the predicted hypothetical proteins generated by prokka...\n";
my $nrDB = "/export/EC1680U/DataBase/Uniprot/Bacteria/Staphylococcus-taxonomy_1279.fasta";
if(! -e "$faa.blastp.tmp"){
	system("blastp -query $faa -db $nrDB -num_threads 32 -max_target_seqs 1 -max_hsps 1 -outfmt '6 qseqid sseqid stitle' -out $faa.blastp.tmp");
}
my $check_blastp = `grep BLASTP $faa.blastp.tmp | wc -l`; chomp($check_blastp);
if($check_blastp == 0){
	system("sed -i '1i Prokka_ID\trefID\tBLASTP' $faa.blastp.tmp");
}
system('perl -pe "s/(sp|tr)\|//" '."$faa.blastp.tmp".' | perl -pe "s/\|(\w+)\t/\t/" > blastp.tmp');
system("sed 's/MULTISPECIES: //g' blastp.tmp > $faa.blastp.tmp");



### Generate report from Prokka annotation ###
print CurrTime() . "Generate report from Prokka annotation...\n";
system("grep -v '##' $prokka_dir/$prefix.gff | grep -P '^TRINITY_DN' | grep -P -v \"\tgene\t\" > $prokka_dir/$prefix.gff.tmp");

# Parse ffn & faa sequences and merge into func_anno.txt
my $oneline_ffn = 'awk \'/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}\' < '.$ffn.' > '.$ffn.'.oneline.ffn';
system("$oneline_ffn");
`sed '/^ *\$/d' $ffn.oneline.ffn | cut -d' ' -f1 |  awk '{printf "%s%s",\$0,(NR%2?FS:RS)}' | sed 's/^>//' | sed 's/ /\t/' > $ffn.oneline.ffn.tmp`;
system("sed -i '1i Prokka_ID\tNucl seq' $ffn.oneline.ffn.tmp");

my $oneline_faa = 'awk \'/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}\' < '.$faa.' > '.$faa.'.oneline.faa';
system("$oneline_faa");
`sed '/^ *\$/d' $faa.oneline.faa | cut -d' ' -f1 |  awk '{printf "%s%s",\$0,(NR%2?FS:RS)}' | sed 's/^>//' | sed 's/ /\t/' > $faa.oneline.faa.tmp`;
system("sed -i '1i Prokka_ID\tProt seq' $faa.oneline.faa.tmp");


open IN, "<$prokka_dir/$prefix.gff.tmp";
open OUT, ">$anno_dir/$prefix\_func_anno.txt.tmp";

while(<IN>){
    chomp;
    my @sp = split/\t/, $_;
    print OUT "$sp[0]\t$sp[1]\t$sp[2]\t$sp[3]\t$sp[4]\t$sp[6]\t";
    my $anno = $sp[8];
    if($anno =~ /ID=(.+);Parent/){
        print OUT "$1\t";
    } else {print OUT "\t";}
    if($anno =~ /eC_number=(.+);Name/){
        print OUT "$1\t";
    } else {print OUT "\t";}
    if($anno =~ /gene=(.+);inference/){
        print OUT "$1\t";
    } else {print OUT "\t";}
    if($anno =~ /product=(.+)$/){
        print OUT "$1\n";
    } else {print OUT "\t\n";}
}
system("sed -i '1i transcript_id\tSource\tFeature\tStart\tEnd\tStrand\tProkka_ID\tEC_number\tGene\tProduct' $anno_dir/$prefix\_func_anno.txt.tmp");


close IN;
close OUT;

my $prokka_table = "$anno_dir/$prefix\_func_anno.txt.tmp";
my $switch_table = "$anno_dir/$prefix\_func_anno.switch.tmp";
my $norepeat_table = "$anno_dir/$prefix\_func_anno.txt.mod";
my $blastp_table = "$faa.blastp.tmp";
my $nucl_seq = "$ffn.oneline.ffn.tmp";
my $prot_seq = "$faa.oneline.faa.tmp";
my $GO_table = "$go_dir/GeneID_GO.txt";
my $COG_table = "$cog_dir/GeneID_COG.txt";
my $report_out = "$anno_dir/$prefix\_func_anno.txt";
my $pathway_table = "$kegg_dir/GeneID_pathway.txt";
my $card_table = "$anno_dir/$prefix\_card.txt";

system("awk -F'\t' 'BEGIN{OFS=FS;}{t=\$1; \$1=\$7; \$7=t;}1' $prokka_table > $switch_table");
system("grep -v repeat_region $switch_table > $norepeat_table");

my $report_cmd = 'prokka <- read.table("'.$norepeat_table.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'blastp <- read.table("'.$blastp_table.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'nucl_seq <- read.table("'.$nucl_seq.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'prot_seq <- read.table("'.$prot_seq.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'cog <- read.table("'.$COG_table.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'go <- read.table("'.$GO_table.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'pathway <- read.table("'.$pathway_table.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'merge1 <- merge(prokka, blastp, by="Prokka_ID", all.x=TRUE);';
$report_cmd .= 'merge2 <- merge(merge1, nucl_seq, by="Prokka_ID", all.x=TRUE);';
$report_cmd .= 'merge3 <- merge(merge2, prot_seq, by="Prokka_ID", all.x=TRUE);';
$report_cmd .= 'merge4 <- merge(merge3, cog, by="Prokka_ID", all.x=TRUE);';
$report_cmd .= 'merge5 <- merge(merge4, go, by="Prokka_ID", all.x=TRUE);';
$report_cmd .= 'merge6 <- merge(merge5, pathway, by="Prokka_ID", all.x=TRUE);';
$report_cmd .= 'write.table(merge6, file="'.$report_out.'", sep="\t", quote=F, na="-", row.names = FALSE);';
system("Rscript" . " -e '" . $report_cmd . "'");

system("$tab2xlsx $report_out $prefix\_func_anno.xlsx");

system("rm -rf $prokka_dir/$prefix.log");
system("rm -rf $prefix\_interpro.* temp");
system("rm -rf $prokka_dir/$prefix.ec");

system("mv 4_Annotation/prokka_func_anno.txt; rm -rf 1_Assembly 2_Subreads 5_Report blastp.tmp Trinity.95.fasta 4_Annotation/prokka_func_anno*; mv 4_Annotation/* .; rm -rf 4_Annotation; mv 3_Prokka Prokka");


print CurrTime() . "Complete!\n";

