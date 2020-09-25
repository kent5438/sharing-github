#! /usr/bin/env perl
### Just only use for BACTERIA ONLY	###
### Need to work by SPAdes assembly ###

use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use Time::Seconds;


my $pwd = `pwd -L`;
chomp($pwd);
my $USER = `id -u -n`;
chomp($USER);
my $rscript = "/usr/bin/Rscript";
my $tab2xls = "/export/EC1680U/perl/bin/tab2xls.pl";
my $tab2xlsx = "/export/EC1680U/perl/bin/tab2xlsx.pl";

my ($help, $contigs, $prefix, $code);
GetOptions(
	"help!"	=> \$help,
	"contigs=s"	=> \$contigs,
	"prefix=s"	=> \$prefix,
	"kegg=s"	=> \$code,
);

if($help){
	die "\n### USAGE ###
* Bacteria amplicon de-novo assembly & annotation pipeline by following tools
- Prokka (prokka)
- COG (rpsblast & cdd2cog.pl)
- GO (interproscan)
- KEGG (ec2kegg.pl)
- Roary (phylogenetic analysis)
- Variant calling (samtools mpileup)

E.g.)
% $0 -p 2009AD47 -c contigs.fasta -k ead

--contigs/-c : Assembled genome by assembler
--prefix/-p : output prefix
--kegg/-k : KEGG_organism_code (http://www.genome.jp/kegg/catalog/org_list.html)
"
} elsif (! $contigs){
	die "\n### ERROR: Please give the contigs/amplicon fasta file\n";
} elsif (! $prefix){
	die "\n### ERROR: Sample name is necessary\n";
} elsif (! $code){
	die "\n### ERROR: KEGG organism code is necessary\n";
} else {
	open OUT, ">>CMD_history.txt";
    print OUT CurrTime() . "$0 -c $contigs -p $prefix -k $code\n";
    close OUT;
}

### Create folders
print CurrTime() . "Create folder in $pwd...\n";
my @folders = qw(1_Assembly 2_Prokka 3_Annotation 4_Report);
my ($assembly_dir, $prokka_dir, $anno_dir, $report_dir) = @folders;
foreach my $folder (@folders){
	if (! -e "$folder"){mkdir "$folder";}
}
my $go_dir = "$anno_dir/GO";
my $kegg_dir = "$anno_dir/KEGG";
my $cog_dir = "$anno_dir/COG";


### Pre-process contig.fasta
print CurrTime() . "Pre-process contig.fasta file...\n";
# Modified bizarre contig header from SPAdes output if necessary
system("perl -i -pe 's/_component_(.+)\$//g' $contigs");

# Contig length distribution
system("/export/EC1680U/perl/bin/Amplicon/Ampliconlength_dist.pl $contigs");

### Generate assembly stats by QUAST-4.5
print CurrTime() . "Generate assembly stats...\n";
my $quast = "python2 /export/EC1680U/software/quast-4.5/quast.py";
my $quast_cmd = "$quast -o $assembly_dir/AssemblyStats $contigs";
if (! -e "$assembly_dir/AssemblyStats"){
	system("$quast_cmd");
}
### Prokka
print CurrTime() . "Processing prokka job...\n";
my $prokka = "/export/EC1680U/software/anaconda2/bin/prokka";
my $prokka_cmd = "$prokka --outdir $prokka_dir --force --prefix $prefix --addgenes --kingdom Bacteria --gcode 11 --cpus 32 --rfam --rnammer $contigs";
if(! -e "$prokka_dir/$prefix.gff"){
	system("$prokka_cmd");
}


### Retrieve eC number from gff, used for the following ec2kegg
print CurrTime() . "Retrieve eC number from gff...\n";
my $ec_cmd = "grep 'eC_number=' $prokka_dir/$prefix.gff | cut -f9 | awk -F';' '{print \$3}' | sed 's/eC_number=//g' > $prokka_dir/$prefix.ec";
system("$ec_cmd");

my $gene_cmd = "grep 'eC_number=' $prokka_dir/$prefix.gff | cut -f9 | awk -F';' '{print \$1}' | sed 's/ID=//g' > $prokka_dir/$prefix.geneID.tmp";
my $paste_cmd = "paste $prokka_dir/$prefix.geneID.tmp $prokka_dir/$prefix.ec > $go_dir/GeneID_ec.txt";
mkdir("$go_dir");
system("$gene_cmd; $paste_cmd");

### GO annotation work by Interproscan
print CurrTime() . "Start GO annotation work...\n";

# Interproscan work (we just want to parse GO data)
print CurrTime() . "Start Interproscan work... (SLOW!)\n";
my $interpro_docker_basal = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` daniel/interproscan:latest";
my $interpro_cmd = "$interpro_docker_basal /interproscan-5.23-62.0/interproscan.sh ".
"-i $prokka_dir/$prefix.faa -appl Pfam -goterms -pa -b $prefix\_interpro -dp";	
if(! -e "$go_dir/GeneID_GO.txt"){
	system("$interpro_cmd");
	system("sudo chown -R $USER:$USER temp $prefix\_interpro.*");
	system("awk -F\"\t\" '{print \$1 \"\t\" \$14}' $prefix\_interpro.tsv > GeneID_GO.txt.tmp");
	my $cmd = 'df <- read.table("GeneID_GO.txt.tmp", header=F, sep="\t");';
	$cmd .= 'df <- df[!(is.na(df$V2) | df$V2==""), ];';
	$cmd .= 'a <- aggregate(V2~V1, data=df, paste, collapse="|");';
	$cmd .= 'write.table(a, "'.$go_dir.'/GeneID_GO.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE);';
	system($rscript . " -e '" . $cmd . "'");
	system("sed -i 's/|/,/g' $go_dir/GeneID_GO.txt");
#	system("sed -i -re 's/_[0-9]+\t/\t/g' $go_dir/GeneID_GO.txt");
	system("sed -i '1i Prokka_ID\tGOs' $go_dir/GeneID_GO.txt");
	system("rm -rf *_interpro.* GeneID_GO.txt.tmp temp");
}


# Parse GO slim (generic) by Trinotate
my $trinotate_home = "/export/EC1680U/software/Trinotate-3.0.2";
my $go_slim_mapping = "$trinotate_home/util/gene_ontology/Trinotate_GO_to_SLIM.pl ".
"$go_dir/GeneID_GO.txt > $go_dir/GO_slim_mapping.txt";
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
system($rscript . " -e '" . $cmd . "'");

# Plot GO barchart
print CurrTime() . "Plot GO mapping graph...\n";
my $GO_mapping_R = "Rscript /export/EC1680U/Rscript/Pacbio/go_mapping_bar_chart.R ".
"$go_dir/GO_mapping.txt $go_dir";
if(! -e "$go_dir/GO_barcahrt.png"){system("$GO_mapping_R");}
system("rm -rf $go_dir/GO_slim_mapping.txt $go_dir/GO_mapping.txt.tmp");

my @go_txt = `ls $go_dir/*.txt`;
chomp(@go_txt);
foreach my $each_go_txt (@go_txt){
    system("$tab2xlsx $each_go_txt $each_go_txt.xlsx");
}


### ec2kegg annotation work...
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
"$kegg_dir/$prefix\_ec2kegg.txt $kegg_dir $code";
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



### COG annotation work...
if(! -e $cog_dir){mkdir("$cog_dir");}
print CurrTime() . "COG annotation work... [rpsbalst stage]\n";
my $cog_db = "/export/EC1680U/DataBase/COG";
my $rpsblast = "/export/EC1680U/software/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/rpsblast";
my $rpsblast_cmd = "$rpsblast -query $prokka_dir/$prefix.faa ".
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
"$cog_dir/func_stats.txt $cog_dir";
if(! -e "$cog_dir/cog_barchart.png"){
    system("$plot_cog_cmd");
}

if (! -e "$cog_dir/GeneID_COG.txt"){
	system("awk -F\"\t\" '{print \$1 \"\t\" \$2}' $cog_dir/protein-id_cog.txt > $cog_dir/protein-id_cog.txt.tmp");
	my $cmd = 'df <- read.table("'.$cog_dir.'/protein-id_cog.txt.tmp", header=F, sep="\t");';
	$cmd .= 'df <- df[!(is.na(df$V2) | df$V2==""), ];';
	$cmd .= 'a <- aggregate(V2~V1, data=df, paste, collapse=",");';
	$cmd .= 'write.table(a, "'.$cog_dir.'/GeneID_COG.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE);';
	system($rscript . " -e '" . $cmd . "'");
	system("sed -i '1i Prokka_ID\tCOGs' $cog_dir/GeneID_COG.txt");
	system("rm -rf $cog_dir/protein-id_cog.txt.tmp");
}

my @cog_txt = `ls $cog_dir/*.txt`;
chomp(@cog_txt);
foreach my $each_cog_txt (@cog_txt){
    system("$tab2xlsx $each_cog_txt $each_cog_txt.xlsx");
}

### Generate contig summary ###
print CurrTime() . "Generate contig summary report...\n";
my $fastalength = "/export/EC1680U/software/exonerate-2.2.0-x86_64/bin/fastalength";
system($fastalength.' '.$contigs.' | sed \'s/ /\t/g\' | awk \'{ t=$1; $1=$2; $2=t; print; }\' | sed \'s/\s/\t/g\' > fastalength.tmp');

# convert a multi-line fasta to a single-line fasta
my $oneline_fa = 'awk \'/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}\' < '.$contigs.' > '.$prefix.'_polished_fixed.oneline.fasta';
system("$oneline_fa");

# split one single-line fasta into multiple single fasta files
my $split_fa = 'cat '.$prefix.'_polished_fixed.oneline.fasta | awk \'/^>/ {if(N>0) printf("\n"); printf("%s\n",$0);++N;next;} { printf("%s",$0);} END {printf("\n");}\' | split -l 2 -a 3 - seq_';
system("$split_fa");

open OUT2, ">GCpct.tmp";
my @seqs = glob("seq_*");
foreach my $seq (@seqs){
    my $GCpct_cmd = `cat $seq | grep -v ">" | awk 'BEGIN{a=0; c=0; g=0; t=0;} {a+=gsub("A",""); c+=gsub("C",""); g+=gsub("G",""); t+=gsub("T","");} END{print (g+c)/(a+g+c+t)*100 }'`;
    my $GCpct = `printf \"%'.2f\"  $GCpct_cmd`;
    print OUT2 "$GCpct\n";
}
close OUT2;

# Estimate plasmid/chromid/chromosome by BLASTN
print CurrTime() . "Estimate plasmid/chromid/chromosome by BLASTN...\n";
if(! -e "$anno_dir/$prefix\_contig_anno.txt"){
	my $blastn = "/export/EC1680U/software/anaconda2/bin/blastn";
    my $blastn_cmd = "$blastn -query $contigs -db /export/EC1680U/DataBase/blastdb/nt -num_threads 32 ".
	"-outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue score stitle qseq staxids sscinames sskingdoms' -max_target_seqs 1 -max_hsps 1 -out $anno_dir/$prefix\_contig_anno.txt.tmp";
    system("$blastn_cmd");
}

open my $blastn_txt, "<$anno_dir/$prefix\_contig_anno.txt.tmp";
open OUT3, ">$anno_dir/$prefix\_contig_anno.txt";
while(<$blastn_txt>){
    chomp;
	my $species;
	my $cov;
	my @sp = split/\t/, $_;
        
	my @cov_sp = split/_cov_/, $sp[0];
	$cov = $cov_sp[1];
    
	foreach(@sp){
        print OUT3 "$_\t";
    }
    print OUT3 "$cov\t";
    
	if($sp[12] =~ /plasmid/){
        print OUT3 "plasmid\n";
    } elsif($sp[12] =~ /chromid/){
        print OUT3 "chromid\n";
    } else {
        print OUT3 "chr\n";
    }
}
close $blastn_txt;
close OUT3;

$paste_cmd = "paste fastalength.tmp GCpct.tmp > $prefix\_contig_summary.txt";
system("$paste_cmd");

system("sed -i '1i Contig\tLength\tGC%' $prefix\_contig_summary.txt");
system("$tab2xlsx $prefix\_contig_summary.txt $assembly_dir/$prefix\_contig_summary.xlsx");
system("rm -rf *.tmp seq_*");

my $check_point = `grep 'NCBI_Prot_ID' $anno_dir/$prefix\_contig_anno.txt | wc -l`; chomp($check_point);
if($check_point == 0){
	system("sed -i '1i $prefix\_Contig\tNCBI_Prot_ID\tpident\tlength\tqlen\tslen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tdescription\tqseq\ttaxID\tSciName\tkingdom\tcoverage\ttype' $anno_dir/$prefix\_contig_anno.txt");
}
system("$tab2xlsx $anno_dir/$prefix\_contig_anno.txt $anno_dir/$prefix\_contig_anno.xlsx");

### Generate report from Prokka annotation
print CurrTime() . "Generate report from Prokka annotation...\n";
system("grep -v '##' $prokka_dir/$prefix.gff | grep -P \"\tID=\" | grep -P -v \"\tgene\t\" > $prokka_dir/$prefix.gff.tmp");

open IN, "<$prokka_dir/$prefix.gff.tmp";
open OUT, ">$anno_dir/$prefix\_func_anno.txt.tmp";

print OUT "Unitig\tSource\tFeature\tStart\tEnd\tStrand\tProkka_ID\tEC_number\tGene\tProduct\n";
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
    if($anno =~ /product=(.+)/){
        print OUT "$1\n";
    } else {print OUT "\t\n";}
}
close IN;
close OUT;

my $prokka_table = "$anno_dir/$prefix\_func_anno.txt.tmp";
my $switch_table = "$anno_dir/$prefix\_func_anno.switch";
my $GO_table = "$anno_dir/GO/GeneID_GO.txt";
my $COG_table = "$anno_dir/COG/GeneID_COG.txt";
my $report_out = "$anno_dir/$prefix\_func_anno.txt";
my $pathway_table = "$kegg_dir/GeneID_pathway.txt";

system("awk -F'\t' 'BEGIN{OFS=FS;}{t=\$1; \$1=\$7; \$7=t;}1' $prokka_table > $switch_table");

my $report_cmd = 'prokka <- read.table("'.$switch_table.'", header=T, sep="\t", check.names=FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'cog <- read.table("'.$COG_table.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'go <- read.table("'.$GO_table.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'pathway <- read.table("'.$pathway_table.'", header=T, sep="\t", check.names = FALSE, fill=TRUE, quote="\"");';
$report_cmd .= 'merge1 <- merge(prokka, cog, by="Prokka_ID", all.x=TRUE);';
$report_cmd .= 'merge2 <- merge(merge1, go, by="Prokka_ID", all.x=TRUE);';
$report_cmd .= 'merge3 <- merge(merge2, pathway, by="Prokka_ID", all.x=TRUE);';
$report_cmd .= 'write.table(merge3, file="'.$report_out.'", sep="\t", quote=F, na="-", row.names = FALSE);';
system($rscript . " -e '" . $report_cmd . "'");

system("$tab2xlsx $report_out $prefix\_func_anno.xlsx");

system("rm -rf $prokka_dir/$prefix.log $prokka_dir/$prefix.sqn");
system("rm -rf $prokka_dir/$prefix.ec $prokka_dir/$prefix.geneID.tmp $prokka_dir/$prefix.gff.tmp $anno_dir/$prefix\_func_anno.txt.tmp $anno_dir/$prefix\_func_anno.switch");

### Report $assembly_dir, $subread_dir, $prokka_dir, $anno_dir, $report_dir
my $pipeline_dir = "/export/EC1680U/Pipeline/PacbioRS-II";
my $help_dir = "/export/EC1680U/kentchen/Template/Help";

system("rsync -avq --exclude='*.fai' $assembly_dir $report_dir");
system("rsync -avq $prokka_dir $report_dir");
system("rsync -avq -m --include='*.xlsx' --include='*.pdf' --include='*.png' --include='*.xls' --include='*/' --exclude='*' $anno_dir $report_dir");
system("rsync -avq $anno_dir/$prefix\_contig_anno.xlsx $report_dir");
system("rsync -avq $prefix\_func_anno.xlsx $report_dir");
system("rsync -avq $help_dir/Amplicon_Help.pdf $report_dir/Help.pdf");
system("rsync -avq $pipeline_dir/README.txt $prokka_dir");

print CurrTime() . "Complete!\n";


sub CurrTime {
    my $t = localtime;
    my $time = "[".$t->hms."]*** ";
    return "$time";
}
