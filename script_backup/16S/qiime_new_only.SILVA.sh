#! /bin/bash

# basic parameter
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
silva_ref_path='/references/SILVA_132_releas'
silva_otu_ref_fasta="/references/SILVA_132_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna"
silva_otu_ref_aligned_fasta="/references/SILVA_132_release/rep_set_aligned/99/99_alignment.fna"
silva_otu_ref_taxonomy="/references/SILVA_132_release/taxonomy/16S_only/99/raw_taxonomy.txt"
silva_ref_tree="/references/SILVA_132_release/trees/99/99_otus.tre"

# env setting
. /usr/local/qiime_software/activate.sh
export PATH=$PATH:/opt/local/software/fastx_toolkit_0.0.13/bin
export AlienTrimmer='/app_tools/AlienTrimmer.jar'
export PATH=/app_tools/mothur/:$PATH

#-- Start --#
START=$(date +%s)


otu_table_97=B.5_OTU_cluster/97/otu_table.otuID.txt
otu_biom_97=B.5_OTU_cluster/97/otu_table.otuID.biom
otu_filter_table_97=B.5_OTU_cluster/97/otu_table.otuID.filter.txt
otu_filter_biom_97=B.5_OTU_cluster/97/otu_table.otuID.filter.biom
otu_gg_table_97=B.5_OTU_cluster/97/otu_table.greengenesID.txt
otu_gg_biom_97=B.5_OTU_cluster/97/otu_table.greengenesID.biom
otu_tax_dedup_table_97=B.5_OTU_cluster/97/otu_table.taxID.filter.dense.dedup.txt
otu_tax_dedup_biom_97=B.5_OTU_cluster/97/otu_table.taxID.filter.dense.dedup.biom

if [ ! -e "$otu_tax_dedup_table_97" ]; then
    echo "### ERROR: Make sure '$otu_tax_dedup_table_97' existed or not!"
    exit 1
fi


#cd $workPath


echo '    B. statistic_analysis'
rm -rf B_statistic_analysis
echo '        plot bar chart'
statistic_sample_mothur.pl $mapping_file B.3_demultiplex_quilty_filter/seqs.fna B.5_OTU_cluster/97/trimmed.unique.opti_mcc.list B.5_OTU_cluster/97/trimmed.unique.opti_mcc.0.03.cons.taxonomy B_statistic_analysis/statistic_sample.txt
plot_StatisticBar.R B_statistic_analysis/statistic_sample.txt B_statistic_analysis
echo '        plot venn chart'
plot_venn.pl $otu_tax_dedup_table_97 $mapping_file B_statistic_analysis/OTU_venn.png "$venn_group"



echo 'C alpha diversity'
echo 'alpha_diversity:metrics chao1,enspie,shannon,simpson_reciprocal,fisher_alpha,goods_coverage,observed_species' > parameter.txt

rm -rf C_alpha_diversity
if [ ! -e "B.5_OTU_cluster/97/otu_table_summary.txt" ]; then
	biom summarize-table -i $otu_tax_dedup_biom_97 -o B.5_OTU_cluster/97/otu_table_summary.txt
fi
maximum_reads=$(grep "Max" B.5_OTU_cluster/97/otu_table_summary.txt | cut -d" " -f 3 | cut -d"." -f 1)

if [ $maximum_reads -ge "30000" ]; then
    multiple_rarefactions.py -i $otu_tax_dedup_biom_97 -m 1 -x $maximum_reads -s 2500 -n 10 -o C_alpha_diversity/rarefaction/
else
    multiple_rarefactions.py -i $otu_tax_dedup_biom_97 -m 1 -x $maximum_reads -s 250 -n 10 -o C_alpha_diversity/rarefaction/
fi

alpha_diversity.py -i C_alpha_diversity/rarefaction/ -o C_alpha_diversity/alpha_div/ -m chao1,enspie,shannon,simpson_reciprocal,fisher_alpha,goods_coverage,observed_species
collate_alpha.py -i C_alpha_diversity/alpha_div/ -o C_alpha_diversity/alpha_div_collated/
make_rarefaction_plots.py -i C_alpha_diversity/alpha_div_collated/ -m $mapping_file -o C_alpha_diversity/alpha_rarefaction_plots

#rm -rf C_alpha_diversity && alpha_rarefaction.py -i $otu_tax_dedup_biom_97 -m $mapping_file -p parameter.txt -t $silva_ref_tree -o C_alpha_diversity -a -O $thread
# /references/gg_13_8_otus/trees/97_otus.tree

echo '    C.4 rank abundance curve'
mkdir -p C.4_rank_abundance && plot_rank_abundance_graph.py -i $otu_tax_dedup_biom_97 -s '*' -n -x -f png -o C.4_rank_abundance/rank_abundance_curve.png
plot_rank_abundance_graph.py -i $otu_tax_dedup_biom_97 -s '*' -x -f png -o C.4_rank_abundance/rank_abundance_curve_L.png

echo '    C. alpha indices table'
rm -rf C_alpha_indices_table && mkdir C_alpha_indices_table
#alpha_diversity.py -i $otu_biom_95 -t $silva_ref_tree -m shannon,PD_whole_tree,chao1,observed_otus,observed_species -o C_alpha_indices_table/95_alpha_indices_table.txt
#alpha_diversity.py -i $otu_tax_dedup_biom_97 -t $silva_ref_tree -m shannon,PD_whole_tree,chao1,observed_otus,observed_species -o C_alpha_indices_table/97_alpha_indices_table.txt
alpha_diversity.py -i $otu_tax_dedup_biom_97 -t $silva_ref_tree -m chao1,enspie,shannon,simpson_reciprocal,fisher_alpha,goods_coverage,observed_species -o C_alpha_indices_table/97_alpha_indices_table.txt
mv C_alpha_indices_table/97_alpha_indices_table.txt C_alpha_indices_table/alpha_indices_table.txt
#merge_indices_table.pl C_alpha_indices_table/95_alpha_indices_table.txt C_alpha_indices_table/97_alpha_indices_table.txt C_alpha_indices_table/alpha_indices_table.txt


echo 'D taxonomic analysis'
echo '    D.1 Krona'
echo '        create taxa summary - L6, L7'
echo 'summarize_taxa:level 1,2,3,4,5,6,7' > parameter.txt # plot all level
rm -rf D.1_taxa_summary && summarize_taxa_through_plots.py -o D.1_taxa_summary -i $otu_tax_dedup_biom_97 -p parameter.txt
echo '        translate Qiime L6.txt to Krona'
rm -rf D.1_translated && translate_Qiime_to_Krnoa.pl D.1_taxa_summary/otu_table.taxID.filter.dense.dedup_L6.txt D.1_translated
echo '        run Krona'
ktImportText -o D.1_Krona.html D.1_translated/*.txt

echo '    D.2 Heatmap'
mkdir -p D.2_heatmap && cd D.2_heatmap
plot_heatmap.pl $mapping_file $workPath/$otu_tax_dedup_table_97 "$group_title"
cd $workPath

echo '    D.3 Phylogenetic tree'
mkdir -p D.3_phylogenetic_tree && cd D.3_phylogenetic_tree
plot_phylogenetic_tree.pl $workPath/$otu_tax_dedup_table_97 $mothur_otu_ref_aligned_fasta
cd $workPath

echo '    D.4 Metastats'
run_Metastats.pl $otu_tax_dedup_table_97 $mapping_file D.4_metastats "$group_title"



echo 'E beta diversity'
echo '    E.1 (Un)weighted unifrac distance PCA/PCoA'
echo 'beta_diversity:metrics  euclidean,unweighted_unifrac,weighted_unifrac' > parameter.txt
rm -rf E.1_PCA_PCoA && beta_diversity_through_plots.py -i $otu_tax_dedup_biom_97 -o E.1_PCA_PCoA -m $mapping_file -t $silva_ref_tree -p parameter.txt
echo '        distance matrix heatmap'
plot_dm_heatmap.R E.1_PCA_PCoA/weighted_unifrac_dm.txt E.1_PCA_PCoA/weighted_dm_heatmap.png
plot_dm_heatmap.R E.1_PCA_PCoA/unweighted_unifrac_dm.txt E.1_PCA_PCoA/unweighted_dm_heatmap.png

echo '    E.2 Clustering analysis of UPGMA'
rm -rf E.2_UPGMA
multiple_rarefactions_even_depth.py -i $otu_tax_dedup_biom_97 -d $sequence_depth -o E.2_UPGMA/rarefaction/
echo '        weighted'
parallel_beta_diversity.py -i E.2_UPGMA/rarefaction/ -o E.2_UPGMA/weighted_unifrac/rare_dm/ -m weighted_unifrac -t $silva_ref_tree -O $thread
rm -rf E.2_UPGMA/weighted_unifrac/rare_dm/BDIV*
upgma_cluster.py -i E.2_UPGMA/weighted_unifrac/rare_dm/ -o E.2_UPGMA/weighted_unifrac/rare_upgma/
consensus_tree.py -i E.2_UPGMA/weighted_unifrac/rare_upgma/ -o E.2_UPGMA/weighted_unifrac/rare_upgma_consensus.tre
tree_compare.py -s E.2_UPGMA/weighted_unifrac/rare_upgma/ -m E.2_UPGMA/weighted_unifrac/rare_upgma_consensus.tre -o E.2_UPGMA/weighted_unifrac/upgma_cmp/
plot_tree.pl E.2_UPGMA/weighted_unifrac/upgma_cmp/master_tree.tre $mapping_file E.2_UPGMA/weighted_unifrac/weighted_UPGMA_tree.png "$group_title" "Weighted UPGMA tree"
echo '        unweighted'
parallel_beta_diversity.py -i E.2_UPGMA/rarefaction/ -o E.2_UPGMA/unweighted_unifrac/rare_dm/ -m unweighted_unifrac -t $silva_ref_tree -O $thread
rm -rf E.2_UPGMA/unweighted_unifrac/rare_dm/BDIV*
upgma_cluster.py -i E.2_UPGMA/unweighted_unifrac/rare_dm/ -o E.2_UPGMA/unweighted_unifrac/rare_upgma/
consensus_tree.py -i E.2_UPGMA/unweighted_unifrac/rare_upgma/ -o E.2_UPGMA/unweighted_unifrac/rare_upgma_consensus.tre
tree_compare.py -s E.2_UPGMA/unweighted_unifrac/rare_upgma/ -m E.2_UPGMA/unweighted_unifrac/rare_upgma_consensus.tre -o E.2_UPGMA/unweighted_unifrac/upgma_cmp/
plot_tree.pl E.2_UPGMA/unweighted_unifrac/upgma_cmp/master_tree.tre $mapping_file E.2_UPGMA/unweighted_unifrac/unweighted_UPGMA_tree.png "$group_title" "Unweighted UPGMA tree"

echo '    E.3 LEfSe'
rm -rf E.3_LEfSe && run_LEfSe.pl D.1_taxa_summary/otu_table.taxID.filter.dense.dedup_L6.txt $mapping_file "$lefse_group" E.3_LEfSe


echo 'F.diversity statistics'
echo '    F.1 NMDS'
rm -rf F.1_NMDS && mkdir F.1_NMDS
nmds.py -i E.1_PCA_PCoA/euclidean_dm.txt -o F.1_NMDS/PCA_nmds.txt && convert_nmds_to_emperor.py F.1_NMDS/PCA_nmds.txt F.1_NMDS/PCA_emperor.txt && make_emperor.py -i F.1_NMDS/PCA_emperor.txt --map_fp $mapping_file -o F.1_NMDS/PCA
nmds.py -i E.1_PCA_PCoA/weighted_unifrac_dm.txt -o F.1_NMDS/PCoA_weighted_nmds.txt && convert_nmds_to_emperor.py F.1_NMDS/PCoA_weighted_nmds.txt F.1_NMDS/PCoA_weighted_emperor.txt && make_emperor.py -i F.1_NMDS/PCoA_weighted_emperor.txt --map_fp $mapping_file -o F.1_NMDS/PCoA_weighted
nmds.py -i E.1_PCA_PCoA/unweighted_unifrac_dm.txt -o F.1_NMDS/PCoA_unweighted_nmds.txt && convert_nmds_to_emperor.py F.1_NMDS/PCoA_unweighted_nmds.txt F.1_NMDS/PCoA_unweighted_emperor.txt && make_emperor.py -i F.1_NMDS/PCoA_unweighted_emperor.txt --map_fp $mapping_file -o F.1_NMDS/PCoA_unweighted

echo '    F.3 Anosim'
mkdir -p F.3_Anosim
if [ $(group_status.pl $mapping_file "ANOSIM" "$group_title") -eq 1 ]
then
    compare_categories.py -i E.1_PCA_PCoA/euclidean_dm.txt -o F.3_Anosim/PCA -m $mapping_file -c "$group_title" --method anosim --num_permutations $num_permutations
    compare_categories.py -i E.1_PCA_PCoA/weighted_unifrac_dm.txt -o F.3_Anosim/PCoA_weighted -m $mapping_file -c "$group_title" --method anosim --num_permutations $num_permutations
    compare_categories.py -i E.1_PCA_PCoA/unweighted_unifrac_dm.txt -o F.3_Anosim/PCoA_unweighted -m $mapping_file -c "$group_title" --method anosim --num_permutations $num_permutations
fi

echo '    F.4 MRPP'
mkdir -p F.4_MRPP
if [ $(group_status.pl $mapping_file "MRPP" "$group_title") -eq 1 ]
then
    compare_categories.py -i E.1_PCA_PCoA/euclidean_dm.txt -o F.4_MRPP/PCA -m $mapping_file -c "$group_title" --method mrpp --num_permutations $num_permutations
    compare_categories.py -i E.1_PCA_PCoA/weighted_unifrac_dm.txt -o F.4_MRPP/PCoA_weighted -m $mapping_file -c "$group_title" --method mrpp --num_permutations $num_permutations
    compare_categories.py -i E.1_PCA_PCoA/unweighted_unifrac_dm.txt -o F.4_MRPP/PCoA_unweighted -m $mapping_file -c "$group_title" --method mrpp --num_permutations $num_permutations
fi

echo '    F.2/F.5 DCA/RDA/CCA'
plot_3A.pl D.1_taxa_summary/otu_table.taxID.filter.dense.dedup_L7.txt $mapping_file 'Species' 'F_3A'

#-- jsyf modified --#

rm -rf D.1_taxa_summary
summarize_taxa.py -i $otu_tax_dedup_biom_97 -o taxa_summary_plots -a -L 1,2,3,4,5,6,7

plot_taxa_summary.py \
        -i taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L1.txt,taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L2.txt,taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L3.txt,taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L4.txt,taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L5.txt,taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L6.txt,taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L7.txt \
        -l Kingdom,Phylum,Class,Order,Family,Genus,Species \
        -c area,bar -o taxa_summary_plots

cd taxa_summary_plots
for i in otu_table.taxID.filter.dense.dedup_*.txt
do
        sed -i '1d' $i
        sed -i 's/^#OTU ID/Taxon/' $i
done
cd ..

sed -i '/Kingdom/i<tr><td class="normal" colspan=2 style="font-size:14"><a href="./taxaPlot_data.xls" target="_blank"> View full table (.xls) </a></td></tr>' taxa_summary_plots/area_charts.html
sed -i '/Kingdom/i<tr><td class="normal" colspan=2 style="font-size:14"><a href="./taxaPlot_data.xls" target="_blank"> View full table (.xls) </a></td></tr>' taxa_summary_plots/bar_charts.html
mv -f taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L* taxa_summary_plots/raw_data/

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

cp B.5_OTU_cluster/97/trimmed.unique.opti_mcc.0.03.cons.taxonomy report_html/data/OTU_Statistic/otu.rep.taxonomy.csv
cp B.5_OTU_cluster/97/trimmed.unique.opti_mcc.0.03.rep.fasta report_html/data/OTU_Statistic/otu.rep.fasta
cp $otu_filter_table_97 report_html/data/OTU_Statistic/
cp $otu_tax_dedup_table_97 report_html/data/OTU_Statistic/

#-- End --#

END=$(date +%s)
DIFF=$(( $END - $START ))
echo ""
echo "Total cost: $DIFF seconds"
