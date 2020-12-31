#! /bin/bash

# basic parameter
ID=`id -u -n`
workPath=$PWD
mapping_file="$workPath/data/Mapping.txt" # absolute path better, files must be fastq format.
group_title="Description"
venn_group="$group_title"
lefse_group="$group_title"

# FLASH v1.2.11 parameter
min_overlap=10       # -m,  The minimum required overlap length between two reads to provide a confident overlap.  Default: 10bp.
max_overlap=65       # -M,  Maximum overlap length expected in approximately 90% of read pairs.  It is by default set to 65bp, which works well for 100bp reads generated from a 180bp library.
mismatch_ratio=0.25  # -x,  Maximum allowed ratio between the number of mismatched base pairs and the overlap length. Default: 0.25.

# AlienTrimmer v0.4 parameter
trim_fasta="$workPath/data/primer.fa" # -c,  for single-end, pair-end (absolute path better, leave '-' for none)
trim_forward_fasta="-"                # -cf, for pair-end only (absolute path better, leave '-' for none)
trim_reverse_fasta="-"                # -cr, for pair-end only (absolute path better, leave '-' for none)
kmer=10;         # -k, k-mer decomposition, between 5 and 15, default 10
mismatches=5;    # -m, maximum mismatches, between 0 and 15, default kmer/2
read_length=150; # -l, minimum read length, default 15
phred_quality=1; # -q, Phred quality score cut-off, between 0 and 40, default 20
percentage=50;   # -p, minimum allowed percentage of correctly called nucleotides, between 0 and 100, default 0

# Mothur parameters
classify_cutoff=60

# other parameters
phred_offset=33
quality=19
num_permutations=999
sequence_depth=110
thread=32

# reference
chimeras_ref='/references/gold.fa'
otu_ref_path='/references/gg_13_8_otus'
otu_ref_tree="$otu_ref_path/trees/97_otus.tree"
mothur_otu_ref_path='/references/gg_13_8_99'
mothur_otu_ref_map="$otu_ref_path/otus/99_otu_map.txt"
mothur_otu_ref_fasta="$mothur_otu_ref_path/gg_13_8_99.fasta"
mothur_otu_ref_aligned_fasta="$mothur_otu_ref_path/gg_13_8_99.refalign"
mothur_otu_ref_taxonomy="$mothur_otu_ref_path/gg_13_8_99.gg.tax"

align_fasta=B.5_OTU_cluster/trimmed.unique.align
otu_table_97=B.5_OTU_cluster/97/otu_table.greengenesID.txt
otu_biom_97=B.5_OTU_cluster/97/otu_table.greengenesID.biom

# env setting
. /usr/local/qiime_software/activate.sh
export PATH=$PATH:/opt/local/software/fastx_toolkit_0.0.13/bin
export AlienTrimmer='/app_tools/AlienTrimmer.jar'
export PATH=/app_tools/mothur/:$PATH



#-- Start --#
START=$(date +%s)


#-- jsyf modified --#

rm -rf D.1_taxa_summary
summarize_taxa.py -i $otu_biom_97 -o taxa_summary_plots -a -L 1,2,3,4,5,6,7

plot_taxa_summary.py \
        -i taxa_summary_plots/otu_table.greengenesID_L1.txt,taxa_summary_plots/otu_table.greengenesID_L2.txt,taxa_summary_plots/otu_table.greengenesID_L3.txt,taxa_summary_plots/otu_table.greengenesID_L4.txt,taxa_summary_plots/otu_table.greengenesID_L5.txt,taxa_summary_plots/otu_table.greengenesID_L6.txt,taxa_summary_plots/otu_table.greengenesID_L7.txt \
        -l Kingdom,Phylum,Class,Order,Family,Genus,Species \
        -c area,bar -o taxa_summary_plots

cd taxa_summary_plots
for i in otu_table.greengenesID_*.txt
do
        sed -i '1d' $i
        sed -i 's/^#OTU ID/Taxon/' $i
done
cd ..

sed -i '/Kingdom/i<tr><td class="normal" colspan=2 style="font-size:14"><a href="./taxaPlot_data.xls" target="_blank"> View full table (.xls) </a></td></tr>' taxa_summary_plots/area_charts.html
sed -i '/Kingdom/i<tr><td class="normal" colspan=2 style="font-size:14"><a href="./taxaPlot_data.xls" target="_blank"> View full table (.xls) </a></td></tr>' taxa_summary_plots/bar_charts.html
mv -f taxa_summary_plots/otu_table.greengenesID_L* taxa_summary_plots/raw_data/

mkdir D.1_taxa_summary

mv taxa_summary_plots D.1_taxa_summary/
cp D.1_taxa_summary/taxa_summary_plots/raw_data/* D.1_taxa_summary/

LINE=`grep -n 'If the lines' C_alpha_diversity/alpha_rarefaction_plots/rarefaction_plots.html | awk -F':' '{print $1}'`
NEXT_LINE=$(($LINE+1))
STRING='<br><br><span style="color: black; font-family:Arial,Verdana; font-size:14; font-weight:bold;"><a href="./rarePlot_data.xls" target="_blank"> View full table (.xls) </a></span>'
sed -i "${NEXT_LINE}i $STRING" C_alpha_diversity/alpha_rarefaction_plots/rarefaction_plots.html

perl /export/EC1680U/perl/bin/16S/rarePlotParser_forCompassPipe.pl C_alpha_diversity/alpha_div_collated/ C_alpha_diversity/alpha_rarefaction_plots/
perl /export/EC1680U/perl/bin/16S/taxaPlotParser_forCompassPipe.pl D.1_taxa_summary/taxa_summary_plots/raw_data/ D.1_taxa_summary/taxa_summary_plots/

#-- jsyf end --#


echo 'Output report html/pdf'
plot_report.rb $quality "report_html" . 'mothur_without_qiime'

cp B.5_OTU_cluster/97/trimmed.unique.an.unique_list.unique.cons.taxonomy report_html/data/OTU_Statistic/otu.rep.taxonomy.csv
cp B.5_OTU_cluster/97/trimmed.unique.an.unique_list.unique.rep.fasta report_html/data/OTU_Statistic/otu.rep.fasta
cp B.5_OTU_cluster/97/otu_table.greengenesID.txt report_html/data/OTU_Statistic/
cp B.5_OTU_cluster/97/otu_table.otuID.txt report_html/data/OTU_Statistic/

#-- End --#
END=$(date +%s)
DIFF=$(( $END - $START ))
echo ""
echo "Total cost: $DIFF seconds"
