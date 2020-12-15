#! /usr/bin/env perl
# Author: kentchen
# Create date: 2017-06-12
# Last update: 2017-10-02

#use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use Time::Seconds;


my $USER = `id -u -n`;
chomp($USER);
GetOptions(
	"input=s"	=> \$input,
	"type=s"	=> \$kingdom,
    "kegg=s"	=> \$code,
	"comp=s"	=> \$comp,
);
if(!($kingdom =~ /euk|bac|arc/)){
	die "\n### ERROR: The organism kingdom need to be euk|bac|arc !\n";
} elsif (! $code){
	die "\n### ERROR: KEGG organism code is necessary\n";
} elsif (! $input){
	die "\n### ERROR: The input file is necessary\n";
} elsif (! $comp){
	die "\n### ERROR: The comparison file is necessary\n";
} else {
    open OUT, ">>CMD_history.txt";
    print OUT CurrTime() . "$0 -i $input -c $comp -t $kingdom -k $code\n";
    close OUT;
}

my $conda_home = "/export/EC1680U/software/anaconda2/bin";
my $perl = "$conda_home/perl";
my $Trinity_home = "/usr/local/src/trinityrnaseq-Trinity-v2.4.0";
my $Trinotate_home = "/export/EC1680U/software/Trinotate-3.0.2";
my $docker_basal_cmd = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` trinityrnaseq/trinityrnaseq:latest";
my $rscript = "/usr/bin/Rscript";
my $tab2xls = "/export/EC1680U/perl/bin/tab2xls.pl"; # need to use for ec2kegg.txt
my $tab2xlsx = "/export/EC1680U/perl/bin/tab2xlsx.pl";
my $qualimap_cmd = "$conda_home/qualimap";

### Create folder
my @folders = qw(1_AssemblyStats 2_Quantitation 3_DiffExpression 4_Annotation);
my ($stat_dir, $quant_dir, $diff_dir, $anno_dir) = @folders;
my $go_dir = "$anno_dir/GO";
my $kegg_dir = "$anno_dir/KEGG";
my $cog_dir = "$anno_dir/COG";

mkdir("$anno_dir");

### Always renew the Trinotate.sqlite (avoid the bug I cannot resolve)
print CurrTime() . "Copy Trinotate.sqlite to your working directory...\n";
system("sudo rsync -av /export/EC1680U/DataBase/Trinotate/Trinotate.sqlite .");
system("sudo chown $USER:$USER Trinotate.sqlite");

### CDS prediction by TransDecoder
print CurrTime() . "Start predicting ORF sequence...\n";
`export PERL5LIB=""`; # cannot find the bug that why PERL5LIB module of DB_file.pm will cause Transdecoder crash, thus empty PERL5LIB variable.
if(! -e "$anno_dir/transdecoder"){mkdir "$anno_dir/transdecoder";}
my $transdec_home = "/export/EC1680U/software/TransDecoder-3.0.1";
my $transdec_long_cmd = "$transdec_home/TransDecoder.LongOrfs -t $stat_dir/Trinity.95.fasta";
my $transdec_predict_cmd = "$transdec_home/TransDecoder.Predict -t $stat_dir/Trinity.95.fasta --cpu 32 ";
if(! -e "$anno_dir/transdecoder/Trinity.95.fasta.transdecoder.pep"){
	system("$transdec_long_cmd");
	system("$transdec_predict_cmd");
	system("mv *.transdecoder* $anno_dir/transdecoder");
}
my $transdec_pep = `ls $anno_dir/transdecoder/*.transdecoder.pep`;
chomp($transdec_pep);


# Generate Qualimap alignment result
print CurrTime() . "Generate Qualimap alignment result...\n";
my @samples = `cat $input | cut -f2`;
chomp(@samples);
foreach $sample (@samples){
	my $gffread_cmd = "/export/EC1680U/software/cufflinks-2.2.1.Linux_x86_64/gffread";
	my $sort_cmd = "/export/EC1680U/software/samtools/samtools sort";
	my $qualimap_bamqc_cmd = "$conda_home/qualimap bamqc";
	if(! -e "$anno_dir/transdecoder/Trinity.95.fasta.transdecoder.gtf"){
		system("$gffread_cmd $anno_dir/transdecoder/Trinity.95.fasta.transdecoder.gff3 -T -o $anno_dir/transdecoder/Trinity.95.fasta.transdecoder.gtf");
	}

	if(! -e "$quant_dir/$sample/bowtie2.sorted.bam"){
		system("$sort_cmd -@ 32 $quant_dir/$sample/bowtie2.bam > $quant_dir/$sample/bowtie2.sorted.bam");
	}

	if(! -e "$quant_dir/$sample/bowtie2.sorted_stats"){
		system("$qualimap_bamqc_cmd -bam $quant_dir/$sample/bowtie2.sorted.bam -gff $anno_dir/transdecoder/Trinity.95.fasta.transdecoder.gtf --java-mem-size=32G");
		my ($line1, $line2, $line3, $line4);
		$line1 = `grep -n '<h2>Input data and parameters' $quant_dir/$sample/bowtie2.sorted_stats/qualimapReport.html | awk -F':' '{print \$1}'`;
		chomp($line1);
		$line2 = `grep -n '<h2>Summary<a' $quant_dir/$sample/bowtie2.sorted_stats/qualimapReport.html | awk -F':' '{print \$1}'`;
		chomp($line2);
		$line2 = $line2 - 1;
		$line3 = `grep -n '<h3>Chromosome stats (inside of regions)' $quant_dir/$sample/bowtie2.sorted_stats/qualimapReport.html | awk -F':' '{print \$1}'`;
		chomp($line3);
		$line4 = `grep -n '<h2>Coverage across reference<a' $quant_dir/$sample/bowtie2.sorted_stats/qualimapReport.html | awk -F':' '{print \$1}'`;
		chomp($line4);
		$line4 = $line4 - 10;
		my $replace1 = 'sed -i "'.$line1.','.$line2.'d" '.$quant_dir.'/'.$sample.'/bowtie2.sorted_stats/qualimapReport.html';
		my $replace2 = 'sed -i "'.$line3.','.$line4.'d" '.$quant_dir.'/'.$sample.'/bowtie2.sorted_stats/qualimapReport.html';
		system("$replace2");
		system("$replace1");
	}
}

# Generate QualiMap multiqc report
if(! -e "$quant_dir/multiqc_qualimap.html"){
	my $multiqc = "/export/EC1680U/software/anaconda2/bin/multiqc";
	my $multiqc_cmd = "$multiqc -f ".
	"-c /export/EC1680U/Pipeline/Trinity/multiqc_qualimap.yaml ".
	"-d $quant_dir/*/bowtie2.sorted_stats ".
	"-o $quant_dir";
	system("$multiqc_cmd");
}


### Sequence homology search by blastx/blastp/hmmscan
print CurrTime() . "Sequence homology search by blastx & blastp...\n";
my $uniDB = "/export/EC1680U/DataBase/Trinotate/uniprot_sprot.pep";
my $nrDB = "/export/EC1680U/DataBase/blastdb/nr";
my $gilist = "/export/EC1680U/DataBase/gilist/S_cerevisiae.prot.180504.gilist";

if(! -e "$anno_dir/blast_result"){mkdir "$anno_dir/blast_result";}
my $blastx_cmd = "/export/EC1680U/software/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastx ".
"-query $stat_dir/Trinity.95.fasta ".
"-db $uniDB ".
#"-db $nrDB -gilist $gilist ".
"-num_threads 32 -max_target_seqs 1 -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle' -evalue 1e-6 > $anno_dir/blast_result/uniprot.blastx.outfmt6";
my $blastp_cmd = "/export/EC1680U/software/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastp ".
"-query $anno_dir/transdecoder/Trinity.95.fasta.transdecoder.pep ".
"-db $uniDB ".
#"-db $nrDB -gilist $gilist ".
"-num_threads 32 -max_target_seqs 1 -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle' -evalue 1e-3 > $anno_dir/blast_result/uniprot.blastp.outfmt6";
if(! -e "$anno_dir/blast_result/uniprot.blastx.outfmt6"){system("$blastx_cmd");}
if(! -e "$anno_dir/blast_result/uniprot.blastp.outfmt6"){system("$blastp_cmd");}

### Multiple features of predicted sequence

# hmmer
print CurrTime() . "Start multiple features prediction... [hmmscan]\n";
if(! -e "$anno_dir/features"){mkdir "$anno_dir/features";}
my $hmmer_cmd = "/export/EC1680U/software/hmmer-3.1b2-linux-intel-x86_64/src/hmmscan ".
"--cpu 32 --domtblout $anno_dir/features/PFAM.out -E 0.01 ".
"/export/EC1680U/DataBase/Trinotate/Pfam-A.hmm $transdec_pep 1> /dev/null";
if(! -e "$anno_dir/features/PFAM.out"){system("$hmmer_cmd");}

# signalp
print CurrTime() . "Start multiple features prediction... [signalp]\n";
my $signalp_cmd = "/export/EC1680U/software/signalp-4.1/signalp ".
"-f short -n $anno_dir/features/signalp.out $transdec_pep 1> /dev/null";
if(! -e "$anno_dir/features/signalp.out"){system("$signalp_cmd");}

# tmhmm
print CurrTime() . "Start multiple features prediction... [tmhmm]\n";
my $tmhmm_cmd = '/export/EC1680U/software/tmhmm-2.0c/bin/tmhmm '.
'--short < '.$transdec_pep.' > '.$anno_dir.'/features/tmhmm.out';
if(! -e "$anno_dir/features/tmhmm.out"){system("$tmhmm_cmd; rm -rf TMHMM_*");}

# rnammer
my $rnammer_cmd = "$Trinotate_home/util/rnammer_support/RnammerTranscriptome.pl ".
"--transcriptome $stat_dir/Trinity.95.fasta --org_type $kingdom ".
"--path_to_rnammer /export/EC1680U/software/Rnammer/rnammer";
if(! -e "$anno_dir/features/Trinity.95.fasta.rnammer.gff"){
	system("$rnammer_cmd");
	system("sudo mv Trinity.95.fasta.rnammer.gff $anno_dir/features");
	system("rm -rf tmp.superscaff* transcriptSuperScaffold*");
}

### Loading database and Generating a Trinotate Annotation Report
print CurrTime() . "Loading the database...\n";
my $gene_trans_map_cmd = "$docker_basal_cmd $Trinity_home/util/support_scripts/get_Trinity_gene_to_trans_map.pl ".
"$stat_dir/Trinity.95.fasta > $stat_dir/Trinity.95.fasta.gene_trans_map";
if(! -e "$stat_dir/Trinity.95.fasta.gene_trans_map"){system("$gene_trans_map_cmd");}

my $loading_db = "$Trinotate_home/Trinotate Trinotate.sqlite init ".
"--gene_trans_map $stat_dir/Trinity.95.fasta.gene_trans_map ".
"--transcript_fasta $stat_dir/Trinity.95.fasta ".
"--transdecoder_pep $anno_dir/transdecoder/Trinity.95.fasta.transdecoder.pep";
my $blastx2db = "$Trinotate_home/Trinotate Trinotate.sqlite LOAD_swissprot_blastx $anno_dir/blast_result/uniprot.blastx.outfmt6";
my $blastp2db = "$Trinotate_home/Trinotate Trinotate.sqlite LOAD_swissprot_blastp $anno_dir/blast_result/uniprot.blastp.outfmt6";
my $pfam2db = "$Trinotate_home/Trinotate Trinotate.sqlite LOAD_pfam $anno_dir/features/PFAM.out";
my $signalp2db = "$Trinotate_home/Trinotate Trinotate.sqlite LOAD_signalp $anno_dir/features/signalp.out";
my $tmhmm2db = "$Trinotate_home/Trinotate Trinotate.sqlite LOAD_tmhmm $anno_dir/features/tmhmm.out";
my $rnammer2db = "$Trinotate_home/Trinotate Trinotate.sqlite LOAD_rnammer $anno_dir/features/Trinity.95.fasta.rnammer.gff";
if(! -e "$anno_dir/Trinotate_report.xls"){system("$loading_db; $blastx2db; $blastp2db; $pfam2db; $signalp2db; $tmhmm2db; $rnammer2db");}

my $trinotate_report = "$Trinotate_home/Trinotate Trinotate.sqlite report > $anno_dir/Trinotate_report.xls";
if(! -e "$anno_dir/Trinotate_report.xls"){system("$trinotate_report");}

my $import_name = "$Trinotate_home/util/annotation_importer/import_transcript_names.pl ".
"Trinotate.sqlite $anno_dir/Trinotate_report.xls";
if(! -e "$anno_dir/Trinotate_report.xls"){system("$import_name");}

### GO data retrieval & plot the graph
print CurrTime() . "Start GO annotation work...\n";
if(! -e "$go_dir"){mkdir "$go_dir";}
my $extract_go_terms = "$Trinotate_home/util/extract_GO_assignments_from_Trinotate_xls.pl ".
"--Trinotate_xls $anno_dir/Trinotate_report.xls -T -I > $go_dir/Trinotate_report.xls.gene_ontology";
if(! -e "$go_dir/Trinotate_report.xls.gene_ontology"){system("$extract_go_terms");}
my $go_slim_mapping = "$Trinotate_home/util/gene_ontology/Trinotate_GO_to_SLIM.pl ".
"$go_dir/Trinotate_report.xls.gene_ontology > $go_dir/GO_slim_mapping.txt";
if(! -e "$go_dir/GO_slim_mapping.txt"){system("$go_slim_mapping");}

open IN, "<$go_dir/GO_slim_mapping.txt";
open OUT, ">$go_dir/GO_mapping.txt";
while(<IN>){
    chomp;
    my @sp=split/\t/, $_;
	no warnings 'uninitialized';
    if(($sp[2] ne "biological_process") && ($sp[2] ne "cellular_component") && ($sp[2] ne "molecular_function")){
        print OUT "$sp[0]\t$sp[1]\t$sp[2]\t$sp[3]\t$sp[4]\n";
    }
}
close IN;
close OUT;

system("sed '1i transcript_id\tGOs' $go_dir/Trinotate_report.xls.gene_ontology > $go_dir/Trans_GO.txt");
system("sed -i '1i Class\tGO_ID\tFunction\tCount\tDescription' $go_dir/GO_mapping.txt");

system("$tab2xlsx $go_dir/Trans_GO.txt $go_dir/Trans_GO.xlsx");
system("$tab2xlsx $go_dir/GO_mapping.txt $go_dir/GO_mapping.xlsx");

# GO enrichment
open IN, "<$comp";
print CurrTime() . "Start GO enrichment work...\n";
my $cut_length_cmd = "cat $quant_dir/final_RSEM_isoforms_report.txt | cut -f 1,2 > $quant_dir/RSEM_isoforms_length.txt.tmp";
system("$cut_length_cmd");

my $enrich_count = `ls $go_dir | grep enrich | wc -l`;
my $comp_count = `cat $comp | wc -l`;
chomp($comp_count);

my ($comparisons, $cmd, @reports, @sig_reports);
if($enrich_count != $comp_count){
	while(<IN>){
		chomp;
		my @tab = split/\t/, $_;
	    my $ctrl_count = scalar (split/,/, $tab[0]);
	    my $case_count = scalar (split/,/, $tab[1]);
	    my $ctrls = join '-', (split/,/, $tab[0]);
	    my $cases = join '-', (split/,/, $tab[1]);
	
		my $comparisons = "$tab[0]-$tab[1]";

		mkdir("$go_dir/$comparisons\_enrich");
		my $docker_mod_trans_cmd = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd`/$diff_dir/$comparisons trinityrnaseq/trinityrnaseq:latest";
	    my $analyze_diff_trans_cmd = "$docker_mod_trans_cmd $Trinity_home/Analysis/DifferentialExpression/analyze_diff_expr.pl ".
	    "--matrix Trans.TMM.counts.matrix ".
	    "-P 0.05 -C 1 --examine_GO_enrichment --GO_annots ../../$go_dir/Trinotate_report.xls.gene_ontology ".
	    "--gene_lengths ../../$quant_dir/RSEM_isoforms_length.txt.tmp --include_GOplot";
	
		system("$analyze_diff_trans_cmd");
		system("sudo chown -R $USER:$USER $diff_dir/$comparisons");
		system("mkdir $diff_dir/$comparisons/tmp");
		system("mv $diff_dir/$comparisons/*GOseq* $diff_dir/$comparisons/Rplots.pdf $go_dir/$comparisons\_enrich");
	    system("mv $diff_dir/$comparisons/*.pdf $diff_dir/$comparisons/*.xlsx $diff_dir/$comparisons/Trans.*.counts.matrix $diff_dir/$comparisons/*.list ".
    	"$diff_dir/$comparisons/*.DE_results $diff_dir/$comparisons/*.subset $diff_dir/$comparisons/*edgeR.count_matrix ".
	    "$diff_dir/$comparisons/tmp/");
	    system("ls $diff_dir/$comparisons/* | grep -v tmp | xargs rm -f");
	    system("mv $diff_dir/$comparisons/tmp/* $diff_dir/$comparisons ; rm -rf $diff_dir/$comparisons/tmp");
		system("rm -rf $go_dir/$comparisons\_enrich/__runGOseq.R");
		system("for FILE in $go_dir/$comparisons\_enrich/Trans.raw.counts.matrix.*; do mv \$FILE \$(echo \$FILE | sed s/Trans\.raw\.counts\.matrix\.//); done");

		my @enriches = `ls $go_dir/$comparisons\_enrich/*.enriched`;
		chomp(@enriches);
		foreach $enrich (@enriches){
			system("$tab2xlsx $enrich $enrich.xlsx");
		}
	}
}


# GO plot
print CurrTime() . "Plot GO mapping graph...\n";
my $GO_mapping_R = "Rscript /export/EC1680U/Rscript/Trinity/go_mapping_bar_chart.R ".
"$go_dir/GO_mapping.txt $go_dir";
if(! -e "$go_dir/GO_barcahrt.png"){system("$GO_mapping_R");}
system("rm -rf $go_dir/GO_slim_mapping.txt");

### Convert GO to EC and do ec2kegg program
print CurrTime() . "Start KEGG annotation work...\n";
print CurrTime() . "Convert GO to EC... (Powered by jsyf)\n";
if(! -e "$kegg_dir"){mkdir("$kegg_dir");}
if((! -s "$go_dir/Trans_EC.txt") or (! -e "$kegg_dir/ec_list.txt")){
	my %go2ec;
	open(IN,'</export/EC1680U/DataBase/ec2go/ec2go_170609.txt');
	while (<IN>) {
	  chomp;
	  if ($_ eq '') {
	    next;
	  }
	  my @array = split("\t",$_);
	  push(@{$go2ec{$array[2]}},$array[0]);
	}
	close(IN);
	
	my %ecGot;
	open(OUT,">$go_dir/Trans_EC.txt");
	open(IN,"<$go_dir/Trinotate_report.xls.gene_ontology");
	while (<IN>) {
	  chomp;
	  if ($_ eq '') {
	    next;
	  }
	  my ($gid,$GO_str) = split("\t",$_);
	  my @GOs = split(',',$GO_str);
	  my @ec;
	  foreach my $GO (@GOs) {
	    if (!defined($go2ec{$GO})) { next; }

	    push(@ec,@{$go2ec{$GO}});
	  }
	  @ec = &uniq(@ec);
	  map $ecGot{$_}=1,@ec;
	  print OUT join("\t",$gid,join(',',@ec))."\n";
	}
	close(IN);
	close(OUT);

	open(OUT,">$kegg_dir/ec_list.txt");
	print OUT join("\n",(sort{$a cmp $b} (keys %ecGot)))."\n";
	close(OUT);
}


my $ec2kegg_bin = "/export/EC1680U/software/ec2kegg";
my $ec2kegg_update_cmd = "$ec2kegg_bin/get_definitions.pl";
my $ec2kegg_cmd = "$ec2kegg_bin/ec2kegg.pl $code $kegg_dir/ec_list.txt > $kegg_dir/ec2kegg.txt";
system("$ec2kegg_update_cmd");


# Parse kegg categories to Trans_KEGG.txt
system('awk \'$2!=""\' '.$go_dir.'/Trans_EC.txt > '.$go_dir.'/Trans_EC.noEmpty.txt');

open my $kegg_table, ">$kegg_dir/Trans_KEGG.txt.tmp.tmp";
open my $ec_txt, "<$go_dir/Trans_EC.noEmpty.txt";

my %hash;
while(<$ec_txt>){
    chomp;
    my @sp = split/\t/, $_;
    my ($trans, $ec) = @sp;
    my @ec_list = split/,/, $ec;
    foreach my $an_ec (@ec_list){
        open IN, "<$kegg_dir/ec2kegg.txt";
        $hash{$trans} = $an_ec;
        while(my $line = <IN>){
            chomp($line);
            my @spp = split/\t/, $line;
            if($spp[12] =~ /$hash{$trans}/){
                print $kegg_table "$trans\t$hash{$trans}\t$spp[1]\n";
            }
        }
		close IN;
    }
}

close $ec_txt;
close $kegg_table;

my $kegg_aggr_cmd = 'df1 <- read.table("'.$kegg_dir.'/Trans_KEGG.txt.tmp.tmp", header=F, sep="\t", check.names=F);';
$kegg_aggr_cmd .= 'x <- aggregate(V3~V1, data=df1, paste, collapse=",");';
$kegg_aggr_cmd .= 'write.table(x, file="'.$kegg_dir.'/Trans_KEGG.txt.tmp", sep="\t", row.names=F, quote=F, na="-", col.names=F)';
system($rscript . " -e '" . $kegg_aggr_cmd . "'");


open IN, "<$kegg_dir/Trans_KEGG.txt.tmp";
open OUT, ">$kegg_dir/Trans_KEGG.txt";
while(<IN>){
	chomp;
	my @sp = split/\t/, $_;
	my @spp = split/,/, $sp[1];
    my @unique = uniq(@spp);
    my @join_unique = join(',', @unique);
    print OUT "$sp[0]\t@join_unique\n";
}
close IN;
close OUT;

system("sed '1i transcript_id\tPathway name' $kegg_dir/Trans_KEGG.txt > $kegg_dir/KEGG_report.txt");


# KEGG pathway plot
print CurrTime() . "Parse KEGG annotation from eC number & plot...\n";
system("$ec2kegg_cmd");
system("sed 's/||/\t/g' $kegg_dir/ec2kegg.txt > $kegg_dir/ec2kegg.txt.tmp");
system("grep -v '^00000' $kegg_dir/ec2kegg.txt.tmp > $kegg_dir/ec2kegg.txt");
system("rm -rf $kegg_dir/ec2kegg.txt.tmp");
system("$tab2xls $kegg_dir/ec2kegg.txt $kegg_dir/ec2kegg.xls");

my $plot_kegg_cmd = "Rscript /export/EC1680U/Rscript/Trinity/kegg_bar_chart.R ".
"$kegg_dir/ec2kegg.txt $kegg_dir $code";
system("$plot_kegg_cmd");
system("cp /export/EC1680U/kentchen/Doc/ec2kegg_README.txt $kegg_dir");


# KEGG up/down regulation
open my $comp_txt, "<$comp";
my $kegg_count = `ls $kegg_dir | grep '_kegg' | wc -l`;

if($kegg_count != $comp_count){
    while(my $comp_line = <$comp_txt>){
        chomp($comp_line);
        my @tab = split/\t/, $comp_line;
        my $ctrl_count = scalar (split/,/, $tab[0]);
        my $case_count = scalar (split/,/, $tab[1]);
        my $ctrls = join '-', (split/,/, $tab[0]);
        my $cases = join '-', (split/,/, $tab[1]);

        my $comparisons = "$tab[0]-$tab[1]";

		mkdir("$kegg_dir/$comparisons\_kegg");
# UP
        open my $up, "<$diff_dir/$comparisons/$comparisons\_UP.list";
        if (-e "$kegg_dir/$comparisons\_kegg/GO_UP.list"){system("rm -rf $kegg_dir/$comparisons\_kegg/GO_UP.list");}
        while(my $uplist = <$up>){
			chomp($uplist);
            system("grep '$uplist' $go_dir/Trinotate_report.xls.gene_ontology >> $kegg_dir/$comparisons\_kegg/GO_UP.list");
        }
        close $up;
		
		if((! -e "$kegg_dir/$comparisons\_kegg/Trans_EC_UP.txt") or (! -e "$kegg_dir/$comparisons\_kegg/ec_list_UP.txt")){
		
		my %go2ec;
		open(IN,'</export/EC1680U/DataBase/ec2go/ec2go_170609.txt');
		while (<IN>) {
		  chomp;
		  if ($_ eq '') {
			next;
		  }
		  my @array = split("\t",$_);
		  push(@{$go2ec{$array[2]}},$array[0]);
		}
		close(IN);

		my %ecGot;
		open(OUT,">$kegg_dir/$comparisons\_kegg/Trans_EC_UP.txt");
		open(IN,"<$kegg_dir/$comparisons\_kegg/GO_UP.list");
		while (<IN>) {
			  chomp;
			  if ($_ eq '') {
				next;
			  }
			  my ($gid,$GO_str) = split("\t",$_);
			  my @GOs = split(',',$GO_str);
			  my @ec;
			  foreach my $GO (@GOs) {
				if (!defined($go2ec{$GO})) { next; }

				push(@ec,@{$go2ec{$GO}});
			  }
			  @ec = &uniq(@ec);	
			  map $ecGot{$_}=1,@ec;
			  print OUT join("\t",$gid,join(',',@ec))."\n";
			}
			close(IN);
			close(OUT);

			open(OUT,">$kegg_dir/$comparisons\_kegg/ec_list_UP.txt");
			print OUT join("\n",(sort{$a cmp $b} (keys %ecGot)))."\n";
			close(OUT);
		}
		
		
		$ec2kegg_cmd = "$ec2kegg_bin/ec2kegg.pl $code $kegg_dir/$comparisons\_kegg/ec_list_UP.txt > $kegg_dir/$comparisons\_kegg/ec2kegg_UP.txt";

		print CurrTime() . "Parse KEGG annotation from eC number & plot...\n";
		system("$ec2kegg_cmd");
		system("sed 's/||/\t/g' $kegg_dir/$comparisons\_kegg/ec2kegg_UP.txt > $kegg_dir/$comparisons\_kegg/ec2kegg.txt.tmp");
		system("grep -v '^00000' $kegg_dir/$comparisons\_kegg/ec2kegg.txt.tmp > $kegg_dir/$comparisons\_kegg/ec2kegg_UP.txt");
		system("rm -rf $kegg_dir/$comparisons\_kegg/ec2kegg.txt.tmp");
		system("$tab2xls $kegg_dir/$comparisons\_kegg/ec2kegg_UP.txt $kegg_dir/$comparisons\_kegg/ec2kegg_UP.xls");

		my $plot_kegg_cmd = "Rscript /export/EC1680U/Rscript/Trinity/kegg_bar_chart.R ".
		"$kegg_dir/$comparisons\_kegg/ec2kegg_UP.txt $kegg_dir/$comparisons\_kegg $code";
		system("$plot_kegg_cmd");
		system("mv $kegg_dir/$comparisons\_kegg/kegg_barchart.png $kegg_dir/$comparisons\_kegg/kegg_barchart_UP.png");
		system("$tab2xlsx $kegg_dir/$comparisons\_kegg/Trans_EC_UP.txt $kegg_dir/$comparisons\_kegg/Trans_EC_UP.xlsx");

# DOWN
		open my $down, "<$diff_dir/$comparisons/$comparisons\_DOWN.list";
		
		if (-e "$kegg_dir/$comparisons\_kegg/GO_DOWN.list"){system("rm -rf $kegg_dir/$comparisons\_kegg/GO_DOWN.list");}
		while(my $dnlist = <$down>){
			chomp($dnlist);
			system("grep '$dnlist' $go_dir/Trinotate_report.xls.gene_ontology >> $kegg_dir/$comparisons\_kegg/GO_DOWN.list");
		}
		close $down;
			
		if((! -e "$kegg_dir/$comparisons\_kegg/Trans_EC_DOWN.txt") or (! -e "$kegg_dir/$comparisons\_kegg/ec_list_DOWN.txt")){
		
		my %go2ec;
		open(IN,'</export/EC1680U/DataBase/ec2go/ec2go_170609.txt');
		while (<IN>) {
		  chomp;
		  if ($_ eq '') {
			next;
		  }
		  my @array = split("\t",$_);
		  push(@{$go2ec{$array[2]}},$array[0]);
		}
		close(IN);

		my %ecGot;
		open(OUT,">$kegg_dir/$comparisons\_kegg/Trans_EC_DOWN.txt");
		open(IN,"<$kegg_dir/$comparisons\_kegg/GO_DOWN.list");
		while (<IN>) {
		  chomp;
		  if ($_ eq '') {
			next;
		  }
		  my ($gid,$GO_str) = split("\t",$_);
		  my @GOs = split(',',$GO_str);
		  my @ec;
		  foreach my $GO (@GOs) {
			if (!defined($go2ec{$GO})) { next; }

			push(@ec,@{$go2ec{$GO}});
		  }
		  @ec = &uniq(@ec);	
		  map $ecGot{$_}=1,@ec;
		  print OUT join("\t",$gid,join(',',@ec))."\n";
		}
		close(IN);
		close(OUT);

		open(OUT,">$kegg_dir/$comparisons\_kegg/ec_list_DOWN.txt");
		print OUT join("\n",(sort{$a cmp $b} (keys %ecGot)))."\n";
		close(OUT);
		}
		
		$ec2kegg_cmd = "$ec2kegg_bin/ec2kegg.pl $code $kegg_dir/$comparisons\_kegg/ec_list_DOWN.txt > $kegg_dir/$comparisons\_kegg/ec2kegg_DOWN.txt";

		system("$ec2kegg_cmd");
		system("sed 's/||/\t/g' $kegg_dir/$comparisons\_kegg/ec2kegg_DOWN.txt > $kegg_dir/$comparisons\_kegg/ec2kegg.txt.tmp");
		system("grep -v '^00000' $kegg_dir/$comparisons\_kegg/ec2kegg.txt.tmp > $kegg_dir/$comparisons\_kegg/ec2kegg_DOWN.txt");
		system("rm -rf $kegg_dir/$comparisons\_kegg/ec2kegg.txt.tmp");
		system("$tab2xls $kegg_dir/$comparisons\_kegg/ec2kegg_DOWN.txt $kegg_dir/$comparisons\_kegg/ec2kegg_DOWN.xls");
		system("$tab2xlsx $kegg_dir/$comparisons\_kegg/Trans_EC_DOWN.txt $kegg_dir/$comparisons\_kegg/Trans_EC_DOWN.xlsx");

		$plot_kegg_cmd = "Rscript /export/EC1680U/Rscript/Trinity/kegg_bar_chart.R ".
		"$kegg_dir/$comparisons\_kegg/ec2kegg_DOWN.txt $kegg_dir/$comparisons\_kegg $code";
		system("$plot_kegg_cmd");
		system("mv $kegg_dir/$comparisons\_kegg/kegg_barchart.png $kegg_dir/$comparisons\_kegg/kegg_barchart_DOWN.png");
	}
}
close $comp_txt;

### COG annotation by CDD instead of eggNOG
print CurrTime() . "COG annotation work...[rpsblast to CDD database]\n";
if(! -e $cog_dir){mkdir("$cog_dir");}
my $cog_db = "/export/EC1680U/DataBase/COG";
my $rpsblast = "/export/EC1680U/software/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/rpsblast";
my $rpsblast_cmd = "$rpsblast -query $anno_dir/transdecoder/Trinity.95.fasta.transdecoder.pep ".
"-db $cog_db/Cog -out $cog_dir/rpsblast.out ".
"-evalue 1e-2 -outfmt 6 -num_threads 32";
if(! -e "$cog_dir/rpsblast.out"){system("$rpsblast_cmd");}

print CurrTime() . "COG annotation work... [Parse COG report]\n";
my $cdd2cog = "$cog_db/cdd2cog.pl";
my $cdd2cog_cmd = "perl $cdd2cog -r $cog_dir/rpsblast.out ".
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
my $plot_cog_cmd = "Rscript /export/EC1680U/Rscript/Trinity/cog_bar_chart.R ".
"$cog_dir/func_stats.txt $cog_dir";
system("$plot_cog_cmd");

# Retrieve COG report table
system("grep -v 'COG#' $cog_dir/rps-blast_cog.txt | cut -f1,13,19 > $cog_dir/rps-blast_cog.txt.tmp");
open IN, "<$cog_dir/rps-blast_cog.txt.tmp";
open OUT, ">$cog_dir/cog_report.txt";
print OUT "transcript_id\tCOGs\tCOG_description\n";
while(<IN>){
    chomp;
    my @sp = split/\t/, $_;
    my @spp = split/\.p/, $sp[0];
    print OUT "$spp[0]\t$sp[1]\t$sp[2]\n";
}
close IN;
close OUT;

my @cog_txt = `ls $cog_dir/*.txt`;
chomp(@cog_txt);
foreach my $each_cog_txt (@cog_txt){
	system("$tab2xlsx $each_cog_txt $each_cog_txt.xlsx");
}


### Generate annotation report
print CurrTime() . "Parse [BLAST] report to $anno_dir/blast_result...\n";
# BLAST report
if (! -e "$anno_dir/blast_result/uniprot.blastp.report.xlsx"){
	open IN, "<$anno_dir/blast_result/uniprot.blastp.outfmt6";
	open OUT, ">$anno_dir/blast_result/uniprot.blastp.report.txt";
	print OUT "transcript_id\tUniprotKB_gene\tpident\tlength\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tDescription\n";
	while(<IN>){
	    chomp;
	    my @sp = split/\t/, $_;
	    my @trans_sp = split/\.p/, $sp[0];
	    my $trans_ID = $trans_sp[0];

		print OUT "$trans_ID\t";
	    for(my $i=1; $i<=$#sp; $i++){
	        print OUT "$sp[$i]\t";
	    }
	    print OUT "\n";
	}
	close IN;
	close OUT;
	system("sed -i 's/\t\$//g' $anno_dir/blast_result/uniprot.blastp.report.txt");
#	system("rm -rf $anno_dir/blast_result/uniprot.blastp.outfmt6");
	system("$tab2xlsx $anno_dir/blast_result/uniprot.blastp.report.txt $anno_dir/blast_result/uniprot.blastp.report.xlsx");
}

if (! -e "$anno_dir/blast_result/uniprot.blastx.report.xlsx"){ 
	open IN, "<$anno_dir/blast_result/uniprot.blastx.outfmt6";
	open OUT, ">$anno_dir/blast_result/uniprot.blastx.report.txt";
	print OUT "transcript_id\tUniprotKB_gene\tpident\tlength\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tDescription\n";
	while(<IN>){
	    chomp;
	    print OUT "$_\t";
	    print OUT "\n";
	}
	close IN;
	close OUT;
	system("sed -i 's/\t\$//g' $anno_dir/blast_result/uniprot.blastx.report.txt");
#	system("rm -rf $anno_dir/blast_result/uniprot.blastx.outfmt6");
	system("$tab2xlsx $anno_dir/blast_result/uniprot.blastx.report.txt $anno_dir/blast_result/uniprot.blastx.report.xlsx");
}

# Features report
print CurrTime() . "Parse [PFAM/signalP/TmHMM/COGs/GOs/KOs] report to $anno_dir/ ...\n";
open IN, "<$anno_dir/Trinotate_report.xls";
open OUT, ">$anno_dir/features/features_report.txt";
while(<IN>){
    chomp;
    my @sp = split/\t/, $_;
    print OUT "$sp[1]\t";
# RNAMMER
	if($sp[3] =~ /\`/){
		my @rnammer_sp = split/\`/, $sp[3];
		foreach my $each_rnammer (@rnammer_sp){
			my @rnammers = split/\^/, $each_rnammer;
			print OUT "$rnammers[0]: $rnammers[1] | ";
		}
		print OUT "\t";
	} elsif($sp[3] =~ /\^/) {
		 my @rnammers = split/\^/, $sp[3];
         print OUT "$rnammers[0]: $rnammers[1]\t";
	} else {
		print OUT "$sp[3]\t";
	}
# PFAM
    if($sp[7] =~ /\`/){
        my @pfam_sp = split/\`/, $sp[7];
        foreach my $each_pfam (@pfam_sp){
            my @pfams = split/\^/, $each_pfam;
            print OUT "$pfams[0]: $pfams[2] | ";
        }
        print OUT "\t";
    } elsif($sp[7] =~ /\^/){
        my @pfams = split/\^/, $sp[7];
        print OUT "$pfams[0]: $pfams[2]\t";
    } else {
        print OUT "$sp[7]\t";
    }
# signalp
    if($sp[8] =~ /\^/){
        $sp[8] =~ s/\^/ \/ /g;
        print OUT "$sp[8]\t";
    } else {
        print OUT "$sp[8]\t";
    }
# TmHMM
    if($sp[9] =~ /\^/){
        $sp[9] =~ s/\^/ \/ /g;
        print OUT "$sp[9]\t";
    } else {
        print OUT "$sp[9]\t";
    }
# COG
#    $sp[10] =~ s/eggnog/COGs (eggnog)/g;
#    if($sp[10] =~ /\`/){
#        my @cog_sp = split/\`/, $sp[10];
#        foreach my $each_cog (@cog_sp){
#            my @cogs = split/\^/, $each_cog;
#            print OUT "$cogs[0]: $cogs[1] | ";
#        }
#        print OUT "\t";
#    } elsif($sp[10] =~ /\^/){
#        my @cogs = split/\^/, $sp[10];
#        print OUT "$cogs[0]: $cogs[1]\t";
#    } else {
#        print OUT "$sp[10]\t";
#    }
# GO
#   $sp[12] =~ s/gene_ontology_blast/GOs/g;
#    if($sp[12] =~ /\`/){
#        my @go_sp = split/\`/, $sp[12];
#        foreach my $each_go (@go_sp){
#            my @gos = split/\^/, $each_go;
#            print OUT "$gos[0], $gos[1], $gos[2] | ";
#        }
#        print OUT "\t";
#    } elsif($sp[12] =~ /\^/){
#        my @gos = split/\^/, $sp[12];
#        print OUT "$gos[0], $gos[1], $gos[2]\t";
#    } else {
#        print OUT "$sp[12]\t";
#    }
# KEGG
#    $sp[11] =~ s/Kegg/KOs/g;
#    if($sp[11] =~ /`(KO:K[0-9]+)/){
#        print OUT "$1";
#    } elsif($sp[11] =~ /KOs/){
#        print OUT "$sp[11]";
#    } else {
#        print OUT ".";
#    }
    print OUT "\n";
}
close IN;
close OUT;
system("sed -i 's/\t\$//g' $anno_dir/features/features_report.txt");
system("$tab2xlsx $anno_dir/features/features_report.txt $anno_dir/features/features_report.xlsx");

# EC number (not belongs to Trinotate_report.xls)
print CurrTime() . "Parse [EC number] report to $anno_dir/ ...\n";
open IN, "<$go_dir/Trans_EC.txt";
open OUT, ">$kegg_dir/EC_report.txt";
print OUT "transcript_id\tEC number\n";
while(<IN>){
	print OUT "$_";
}
close IN;
close OUT;

system("$tab2xlsx $kegg_dir/EC_report.txt $kegg_dir/EC_report.xlsx");

system("sudo rm -rf tmp.*");

### Using KO to parse organism-specific pathway (more detail than using ec number)
#print CurrTime() . "Parse organism-specific pathway by KO...\n";

#system("/export/EC1680U/perl/bin/Trinity/ko2pathway.sh $code");



sub CurrTime {
    my $t = localtime;
    my $time = "[".$t->hms."]*** ";
    return "$time";
}

sub uniq {
      my %seen;
      foreach (@_) {
      if ($_ eq '') {next;}

      $seen{$_} = 1;
    }
    my @uniq = sort{$a cmp $b} (keys %seen);
    return(@uniq);
}
