#! /usr/bin/env perl
### Just only use for BACTERIA ONLY	###
### Should be worked after SPAdes assembly ###

use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use Time::Seconds;
use lib '/export/EC1680U/kentchen/lib';
use Time qw(CurrTime);

system(". activate illumina-assembly");

my $pwd = `pwd -L`;
chomp($pwd);
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
* Bacteria illumina de-novo assembly & annotation pipeline by following tools
- shovill (FLASH-SPAdes-Pilon)
- scaffolding & gapfilling (SSPACE-Gapfiller)
- Prokka (prokka)
- COG (rpsblast & cdd2cog.pl)
- GO (interproscan)
- KEGG (ec2kegg.pl)

E.g.) assumed using shovill-SSPACE-gapfiller
% $0 -p 2009AD47 -k ead

--prefix/-p : output prefix
--kegg/-k : KEGG_organism_code (https://www.genome.jp/kegg/catalog/org_list.html)
"
} elsif (! $prefix){
	die "\n### ERROR: Sample name is necessary\n\n";
} elsif (! $code){
	die "\n### ERROR: KEGG organism code is necessary\n\n";
} else {
	open OUT, ">>CMD_history.txt";
    print OUT CurrTime() . "\n$0 -p $prefix -k $code\n\n";
    close OUT;
}

### Create folders
print CurrTime() . "Create folder in $pwd...\n";

#my $check_loc = `echo $pwd | grep '_shovill'`; chomp($check_loc);
#if(! $check_loc){die "\n### ERROR: Please run this script in 'sampleName_shovill' !!!\n\n";}

my @folders = qw(0_ReadQC 1_Assembly 2_Prokka 3_Annotation 4_Report);
my ($QC_dir, $assembly_dir, $prokka_dir, $anno_dir, $report_dir) = @folders;
foreach my $folder (@folders){
	if (! -e "$folder"){mkdir "$folder";}
}
my $go_dir = "$anno_dir/GO";
my $kegg_dir = "$anno_dir/KEGG";
my $cog_dir = "$anno_dir/COG";


### Reads QC results
print CurrTime() . "Running QC result by multiQC...\n";

if (! -e "reads"){
    if (! -e "../reads"){die "### ERROR: Please make sure 'reads/' is in the root folder";}
    mkdir("reads");
    chdir("reads");
    system("ln -s ../../reads/$prefix* .");
    chdir("..");
}

if(! -e "fastqc_result"){
	mkdir("fastqc_result");
	system("fastqc -t 32 reads/$prefix* -o fastqc_result");
}
if(! -e "0_ReadQC/multiqc_report.html"){
#	system("multiqc -c /export/EC1680U/Pipeline/Miseq_Assembly/multiqc_config.yaml -f fastqc_result/ -o 0_ReadQC");
	system("/export/EC1680U/perl/bin/Assembly/fastq-stats-k.pipe.sh 301");
}

### Pre-process contig.fasta
print CurrTime() . "Pre-process contig.fasta file...\n";

# contigs & scaffolds are assumed generating from shovill-sspace-gapfiller
print CurrTime() . "Pre-process contig.fasta file...\n";

# contigs & scaffolds are assumed generating from shovill-sspace-gapfiller
my $contigs = "1_Assembly/contigs.fasta";
my $scaffolds = "1_Assembly/scaffolds.fasta";

if (! -e "$contigs"){
      if(! -e "shovill/contigs.fa"){die "\n### ERROR: Please make sure you have run '% shovill_sspace_gapfiller.sh $prefix'\n\n";}
      else {system("cp shovill/contigs.fa $contigs");}
}
if(! -e "$scaffolds"){
#   if(! -e "standard_output/standard_output.gapfilled.final.fa"){die "\n### ERROR: Please make sure you have run '% shovill_sspace_gapfiller.sh $prefix'\n\n";}
    if(! -e "shovill/standard_output/standard_output.gapfilled.final.fa"){system("ln -s $pwd/shovill/standard_output/standard_output.final.scaffolds.fasta $pwd/shovill/standard_output/standard_output.gapfilled.final.fa");}
    else {system("cp shovill/standard_output/standard_output.gapfilled.final.fa $scaffolds");}
}

# Modified bizarre contig header from SPAdes output if necessary
system("perl -pi -e 's/_component_(.+)\$//g' $contigs");
system("perl -pi -e 's/\\.(.+)\$//g' $contigs");

### Generate assembly stats by QUAST-4.5
print CurrTime() . "Generate assembly stats...\n";
my $quast_cmd = "quast -o $assembly_dir/AssemblyStats $contigs $scaffolds";
if (! -e "$assembly_dir/AssemblyStats"){
	system("$quast_cmd");
}

system("sed -i '1d' $assembly_dir/AssemblyStats/report.tsv");
system("sed -i '1i $prefix\tcontigs\tscaffolds' $assembly_dir/AssemblyStats/report.tsv");
system('sed -i "s/# contigs (/contigs (/g" '.$assembly_dir.'/AssemblyStats/report.tsv');
system('sed -i "s/Total length (/# Total length (/g" '.$assembly_dir.'/AssemblyStats/report.tsv');


### Prokka
print CurrTime() . "Processing prokka job...\n";
my $prokka_cmd = "prokka --outdir $prokka_dir --force --prefix $prefix --addgenes --kingdom Bacteria --gcode 11 --cpus 32 --rfam --rnammer $scaffolds";
if(! -e "$prokka_dir/$prefix.gff"){
	system("$prokka_cmd");
}
system('sed -i "s/organism/features/g" '."$prokka_dir/$prefix.txt");
system('sed -i "s/Genus species strain/count/g" '."$prokka_dir/$prefix.txt");
system('sed -i "s/: /\t/g" '."$prokka_dir/$prefix.txt");


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
my $interpro_docker_basal = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` kentchendocker/interproscan:latest";
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
my $GO_mapping_R = "Rscript /export/EC1680U/Rscript/Pacbio/go_mapping_bar_chart.R ".
"$go_dir/GO_mapping.txt $go_dir &> /dev/null";
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


### COG annotation work...
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

### Generate contig summary ###
print CurrTime() . "Generate contig summary report...\n";
system("fastalength".' '.$scaffolds.' | sed \'s/ /\t/g\' | awk \'{ t=$1; $1=$2; $2=t; print; }\' | sed \'s/ /\t/g\' > fastalength.tmp');

# convert a multi-line fasta to a single-line fasta
my $oneline_fa = 'awk \'/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}\' < '.$scaffolds.' > '.$scaffolds.'.oneline.fasta';
system("$oneline_fa");
system("sed -i '/^ *\$/d' $scaffolds.oneline.fasta"); # remove empty line


# split one single-line fasta into multiple single fasta files
my $split_fa = 'cat '.$scaffolds.'.oneline.fasta | awk \'/^>/ {if(N>0) printf("\n"); printf("%s\n",$0);++N;next;} { printf("%s",$0);} END {printf("\n");}\' | split -l 2 -a 3 - seq_';
system("$split_fa");

open OUT2, ">GCpct.tmp";
my @seqs = glob("seq_*");
foreach my $seq (@seqs){
    my $GCpct_cmd = `cat $seq | grep -v ">" | awk 'BEGIN{a=0; c=0; g=0; t=0;} {a+=gsub("A",""); a+=gsub("a",""); c+=gsub("C",""); c+=gsub("c",""); g+=gsub("G",""); g+=gsub("g",""); t+=gsub("T",""); g+=gsub("t","");} END{print (g+c)/(a+g+c+t)*100 }'`;
    my $GCpct = `printf \"%'.2f\"  $GCpct_cmd`;
    print OUT2 "$GCpct\n";
}
close OUT2;

# Estimate plasmid/chromid/chromosome by BLASTN
print CurrTime() . "Estimate plasmid/chromid/chromosome by BLASTN...\n";
my $nrDB = "/export/EC1680U/DataBase/Uniprot/Bacteria/Staphylococcus-taxonomy_1279.fasta";
#my $gilist = "/export/EC1680U/DataBase/gilist/Archaea.prot.180927.gilist"; # all bacteria gilist
my $ntDB = "/export/EC1680U/DataBase/blastdb/nt";
if(! -e "$anno_dir/$prefix\_contig_anno.txt.tmp"){
    my $blastn_cmd = "blastn -query $scaffolds -db $ntDB -num_threads 32 ".
	"-outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue score stitle staxids sscinames sskingdoms' -max_target_seqs 3 -max_hsps 1 -out $anno_dir/$prefix\_contig_anno.txt.tmp";
    system("$blastn_cmd");
}

my $cov_exist;
open my $blastn_txt, "<$anno_dir/$prefix\_contig_anno.txt.tmp";
open OUT3, ">$anno_dir/$prefix\_contig_anno.txt";
while(<$blastn_txt>){
    chomp;
	my $species;
	my $cov;
    my @sp = split/\t/, $_;
        
    foreach(@sp){
        print OUT3 "$_\t";
    }

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

my $cov_fasta_header = "grep '>' $scaffolds | sed 's/^>//' | awk -F'cov=' '{print \$2}' | cut -d' ' -f1 > cov.tmp";
system("$cov_fasta_header");

$paste_cmd = "paste fastalength.tmp GCpct.tmp cov.tmp > $prefix\_contig_summary.txt";
system("$paste_cmd");

system("join -t '\t' -1 1 $prefix\_contig_summary.txt -2 1 $anno_dir/$prefix\_contig_anno.txt > contig_anno.tmp");
system("cp contig_anno.tmp $anno_dir/$prefix\_contig_anno.txt");

system("sed -i '1i $prefix\_Contig\tLength\tGC%\tCoverage' $prefix\_contig_summary.txt");

system("$tab2xlsx $prefix\_contig_summary.txt $assembly_dir/$prefix\_contig_summary.xlsx");

system("rm -rf *.tmp seq_*");

my $check_point = `grep '_Contig' $anno_dir/$prefix\_contig_anno.txt | wc -l`; chomp($check_point);
if($check_point == 0){
	system("sed -i '1i $prefix\_Contig\tLength\tGC%\tCoverage\trefID\tpident\tlength\tqlen\tslen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tdescription\ttaxID\tSciName\tkingdom\ttype' $anno_dir/$prefix\_contig_anno.txt");
} else {
	system("sed -i '1d' $anno_dir/$prefix\_contig_anno.txt");
	system("sed -i '1i $prefix\_Contig\tLength\tGC%\tCoverage\trefID\tpident\tlength\tqlen\tslen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tdescription\taxID\tSciName\tkingdom\ttype' $anno_dir/$prefix\_contig_anno.txt");
}

system("$tab2xlsx $anno_dir/$prefix\_contig_anno.txt $anno_dir/$prefix\_contig_anno.xlsx");


## Using BLASTP to fill up the predicted hypothetical protein from prokka
print CurrTime() . "Using BLASTP to fill up the predicted hypothetical proteins generated by prokka...\n";
my $faa = "$prokka_dir/$prefix.faa";
my $ffn = "$prokka_dir/$prefix.ffn";
if(! -e "$faa.blastp.tmp"){
#	system("blastp -query $faa -db $nrDB -num_threads 32 -gilist $gilist -max_target_seqs 1 -max_hsps 1 -outfmt '6 qseqid sseqid stitle' -out $faa.blastp.tmp");
	system("blastp -query $faa -db $nrDB -num_threads 32 -max_target_seqs 1 -max_hsps 1 -outfmt '6 qseqid sseqid stitle' -out $faa.blastp.tmp");
}
my $check_blastp = `grep BLASTP $faa.blastp.tmp | wc -l`; chomp($check_blastp);
if($check_blastp == 0){
	system("sed -i '1i Prokka_ID\trefID\tBLASTP' $faa.blastp.tmp");
}
system('perl -pe "s/(sp|tr)\|//" '."$faa.blastp.tmp".' | perl -pe "s/\|(\w+)\t/\t/" > blastp.tmp');
system("sed 's/MULTISPECIES: //g' blastp.tmp > $faa.blastp.tmp");


### Generate report from Prokka annotation
print CurrTime() . "Generate report from Prokka annotation...\n";
system("grep -v '##' $prokka_dir/$prefix.gff | grep -P \"\tID=\" | grep -P -v \"\tgene\t\" | sed 's/;protein_id=*\$//g' > $prokka_dir/$prefix.gff.tmp");

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
        if($anno =~ /product=(.+);protein_id=(.+)$/){
            print OUT "$1\n";
        } elsif($anno =~ /product=(.+)$/) {
            print OUT "$1\n";
        }
    } else {print OUT "\t\n";}
}
system("sed -i '1i Unitig\tSource\tFeature\tStart\tEnd\tStrand\tProkka_ID\tEC_number\tGene\tProduct' $anno_dir/$prefix\_func_anno.txt.tmp");

close IN;
close OUT;

my $prokka_table = "$anno_dir/$prefix\_func_anno.txt.tmp";
my $switch_table = "$anno_dir/$prefix\_func_anno.switch.tmp";
my $norepeat_table = "$anno_dir/$prefix\_func_anno.mod.tmp";
my $blastp_table = "$faa.blastp.tmp";
my $nucl_seq = "$ffn.oneline.ffn.tmp";
my $prot_seq = "$faa.oneline.faa.tmp";
my $GO_table = "$anno_dir/GO/GeneID_GO.txt";
my $COG_table = "$anno_dir/COG/GeneID_COG.txt";
my $report_out = "$anno_dir/$prefix\_func_anno.txt";
my $pathway_table = "$kegg_dir/GeneID_pathway.txt";

system("awk -F'\t' 'BEGIN{OFS=FS;}{t=\$1; \$1=\$7; \$7=t;}1' $prokka_table > $switch_table");
system("grep -v repeat_region $switch_table > $norepeat_table");

my $report_cmd = 'prokka <- read.table("'.$norepeat_table.'", header=T, sep="\t", check.names=FALSE, fill=TRUE, quote="\"");';
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

### Report $assembly_dir, $subread_dir, $prokka_dir, $anno_dir, $report_dir
my $pipeline_dir = "/export/EC1680U/Pipeline/Miseq_Assembly";
my $help_dir = "/export/EC1680U/kentchen/Template/Help";

system("rsync -avq $QC_dir $report_dir");
system("rsync -avq --exclude='*.fai' --exclude='*tmp*' --exclude='*oneline*' $assembly_dir $report_dir");
system("rsync -avq --exclude='*tmp*' --exclude='*oneline*'  $prokka_dir $report_dir");
system("rsync -avq -m --include='*.xlsx' --include='*.pdf' --include='*.png' --include='*.xls' --include='*/' --exclude='*' --exclude='*tmp*' --exclude='*oneline*' $anno_dir $report_dir");
system("rsync -avq $anno_dir/$prefix\_contig_anno.xlsx $report_dir");
system("rsync -avq $prefix\_func_anno.xlsx $report_dir");
system("rsync -avq $help_dir/Assembly_Help.pdf $report_dir/Help.pdf");
system("rsync -avq $pipeline_dir/README.txt $report_dir/$prokka_dir");


### Generate Report by Rmarkdown

# Copy report_template to current dir)
print CurrTime() . "Copy report_template to current dir...\n";
my $report_template = "/export/EC1680U/Pipeline/Miseq_Assembly/report_template";
if(! -e "$report_template"){die "\n### ERROR: Please make sure '$report_template' existed!\n\n";}
system("rsync -avq $report_template .");

# Transfer 5_Report to report_template
print CurrTime() . "Transfer 4_Report/* to report_template/files and running build_site.R...\n";
system("rsync -avq 4_Report/* report_template/files");

# Build R markdown Script
print CurrTime() . "Build Rmarkdown-Style Report...\n";
system("sed -i 's/<prefix>/$prefix/g' report_template/index.Rmd");
system("sed -i 's/<prefix>/$prefix/g' report_template/experiment.Rmd");
system("sed -i 's/<prefix>/$prefix/g' report_template/seqQC.Rmd");
system("sed -i 's/<prefix>/$prefix/g' report_template/assembly.Rmd");
system("sed -i 's/<prefix>/$prefix/g' report_template/annotation.Rmd");
system("sed -i 's/<prefix>/$prefix/g' report_template/tools.Rmd");
system("Rscript /export/EC1680U/Pipeline/Miseq_Assembly/build_site.R");

# Rsync report to root folder & remove report_template
system("rsync -avq report_template/report .");
system("rm -rf report_template blastp.tmp");


print CurrTime() . "Complete!\n";

