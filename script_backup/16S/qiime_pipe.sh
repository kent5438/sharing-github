#! /bin/bash
### deault output: "./qiime_work" ###
### This pipe could be run on sheep only so far ###

output="qiime_work"
username=`id -u -n`

if [ ! -e "$output" ]; then
	mkdir $output
fi
cd $output

cp ../6_biom/*.biom .

qiime="sudo docker run --rm -v `pwd`:`pwd` -w `pwd` bwawrik/qiime:latest"
qiime_taxa="sudo docker run --rm -v `pwd`:`pwd` -w `pwd`/taxa_plots bwawrik/qiime:latest"
rscript="/export/EC1680U/software/anaconda2/bin/Rscript"

### generate the map file (!!! SampleID IS THE FIXED FORMAT for QIIME !!!)
echo "*** Generate the map file..."

if [ -e "otu_table.biom" ]; then
	rm -f otu_table.biom
fi

if [ ! -e "Map_file.txt" ]; then
	echo "#SampleID" > Map_file.txt
	ls *.biom | grep -v "otu_table" | sort -V | sed 's/\.biom//g' >> Map_file.txt
fi

sample_number=$(grep -v "^#SampleID" Map_file.txt | wc -l)


### merge the biom files
echo "*** Merge the biom files..."
#total_biom=$(ls *.biom | grep -v "otu_table" | sort -u -V | sed ':a;N;$!ba;s/\n/,/g')
total_biom=$(sed '1d' Map_file.txt | awk '{print $0".biom"}' | sed ':a;N;$!ba;s/\n/,/g')
$qiime merge_otu_tables.py -i $total_biom -o otu_table.biom
$qiime biom summarize-table -i otu_table.biom -o otu_table_summary.txt
maximum_reads=$(grep "Max" otu_table_summary.txt | cut -d" " -f 3 | cut -d"." -f 1)


### convert biom to otu table & PCA format
echo "*** Convert biom to otu table..."
$qiime biom convert -i otu_table.biom -o otu_table.txt --header-key taxonomy --to-tsv --table-type="OTU table"
cat otu_table.txt | grep -v 'biom' | awk -F'\t' 'BEGIN{OFS="\t"} {$1=$NF;  print}' | awk -F'\t' 'BEGIN{OFS="\t"} {$NF=""; print}' > otu_table.forPCA.txt
sed -i 's/\t$//' otu_table.forPCA.txt
sudo chown -R $username:$username *


### Generate area/bar taxa plot
echo "*** Generate area/bar taxa plot..."
$qiime summarize_taxa.py -i otu_table.biom -o taxa_plots -a  -L 2,3,4,5,6,7,8
$qiime_taxa plot_taxa_summary.py \
	-i otu_table_L2.txt,otu_table_L3.txt,otu_table_L4.txt,otu_table_L5.txt,otu_table_L6.txt,otu_table_L7.txt,otu_table_L8.txt \
	-l Kingdom,Phylum,Class,Order,Family,Genus,Species \
	-c area,bar -o .
sudo chown -R $username:$username *

cd taxa_plots
for i in otu_table*.txt
do
	sed -i '1d' $i
	sed -i 's/^#OTU ID/Taxon/' $i
done
cd ..

sed -i '/Kingdom/i<tr><td class="normal" colspan=2 style="font-size:14"><a href="./taxaPlot_data.xls" target="_blank"> View full table (.xls) </a></td></tr>' taxa_plots/area_charts.html
sed -i '/Kingdom/i<tr><td class="normal" colspan=2 style="font-size:14"><a href="./taxaPlot_data.xls" target="_blank"> View full table (.xls) </a></td></tr>' taxa_plots/bar_charts.html
mv -f taxa_plots/otu_table_L* taxa_plots/raw_data/

ln -s /export/EC1680U/perl/bin/16S/taxaPlotParser.pl .
sudo chown -R $username:$username taxa_plots
perl taxaPlotParser.pl taxa_plots/raw_data taxa_plots


### Generate Rank-abundance curve
echo "*** Rank-abundance curve..."
$qiime plot_rank_abundance_graph.py -i otu_table.biom -s '*' -o rank_abundance_graph.L -x -f pdf
$qiime plot_rank_abundance_graph.py -i otu_table.biom -s '*' -o rank_abundance_graph.L -x -f png
$qiime plot_rank_abundance_graph.py -i otu_table.biom -s '*' -o rank_abundance_graph.L -x -f svg
$qiime plot_rank_abundance_graph.py -i otu_table.biom -s '*' -o rank_abundance_graph -n -x -f pdf
$qiime plot_rank_abundance_graph.py -i otu_table.biom -s '*' -o rank_abundance_graph -n -x -f png
$qiime plot_rank_abundance_graph.py -i otu_table.biom -s '*' -o rank_abundance_graph -n -x -f svg


### Generate heatmap plot
echo "*** Generate genus heatmap plot..."
#$qiime make_otu_heatmap.py -i otu_table.biom -o OTUs_heatmap.pdf (new version of qiime does not support)
grep -v "unclassified" taxa_plots/raw_data/otu_table_L7.txt > heatmap.txt
/export/EC1680U/perl/bin/16S/extract_genus_forHeatmap.pl
sed -i 's/\t$//' heatmap.genus.txt
(head -n 1 heatmap.genus.txt && tail -n +2 heatmap.genus.txt | sort -u -k1,1) > heatmap.genus.sorted.uniq.txt
$rscript /export/EC1680U/Rscript/16S/heatmap_eachSample.R heatmap.genus.sorted.uniq.txt .


### Alpha diversity (rarefaction plot)
echo "alpha diversity (rarefaction plot)..."
if [ $maximum_reads -ge "20000" ]; then
	$qiime multiple_rarefactions.py -i otu_table.biom -m 1 -x $maximum_reads -s 1000 -n 10 -o rare_step_5000div_1
else
	$qiime multiple_rarefactions.py -i otu_table.biom -m 1 -x $maximum_reads -s 100 -n 10 -o rare_step_5000div_1
fi


# ACE seems to not support anymore, remove
$qiime alpha_diversity.py -i rare_step_5000div_1 -o rare_step_5000div_2/ -m chao1,enspie,shannon,simpson_reciprocal,fisher_alpha,goods_coverage,observed_species
$qiime collate_alpha.py -i rare_step_5000div_2/ -o rare_step5000div_3_alpha_collated
$qiime make_rarefaction_plots.py -i rare_step5000div_3_alpha_collated/ -m Map_file.txt -o rare_alpha_plot
sudo chown -R $username:$username *

mkdir meta16S_OUTPUT
mkdir meta16S_OUTPUT/OTUs_heatmap
cp heatmap.genus.sorted.uniq.txt OTUs_heatmap.pdf meta16S_OUTPUT/OTUs_heatmap
mv taxa_plots meta16S_OUTPUT
mv rare_alpha_plot meta16S_OUTPUT
mkdir meta16S_OUTPUT/rank_abundance_graph
rsync -av rank_abundance_graph* meta16S_OUTPUT/rank_abundance_graph

ln -s /export/EC1680U/perl/bin/16S/rarePlotParser.pl .
perl rarePlotParser.pl rare_step5000div_3_alpha_collated/ meta16S_OUTPUT/rare_alpha_plot/

LINE=`grep -n 'If the lines' meta16S_OUTPUT/rare_alpha_plot/rarefaction_plots.html | awk -F':' '{print $1}'`
NEXT_LINE=$(($LINE+1))
STRING='<br><br><span style="color: black; font-family:Arial,Verdana; font-size:14; font-weight:bold;"><a href="./rarePlot_data.xls" target="_blank"> View full table (.xls) </a></span>'
sed -i "${NEXT_LINE}i $STRING" meta16S_OUTPUT/rare_alpha_plot/rarefaction_plots.html

### Beta diversity
echo "beta diversity (correlation / PCA / PCoA-2D)..."
ln -s /export/EC1680U/perl/bin/16S/alphaBetaDiversity.pl .
BETA_DIR=meta16S_OUTPUT/AlphaBetaDiversity
mkdir $BETA_DIR
# need to modified alphaBetaDiversity.pl content for docker distribution
perl alphaBetaDiversity.pl otu_table.biom $BETA_DIR
rm -f 99_otus.tree UNITEv6_sh_99.aln.tre


### Correlation plot
echo "*** Generate correlation heatmap plot..."
$rscript /export/EC1680U/Rscript/16S/correlation_allSample.R $BETA_DIR/weighted_unifrac_otu_table.txt Sample_correlation.pdf
cp Sample_correlation.pdf meta16S_OUTPUT/OTUs_heatmap

# PCA
rsync -avq otu_table.forPCA.txt meta16S_OUTPUT
if [ $(cat Map_file.txt | grep -v '#' | wc -l) -gt 2 ]; then
	$rscript /export/EC1680U/Rscript/16S/PCA_plot.R
fi

# PCoA-2D (plot by weighted unifrac)
#PCoA_DIR="meta16S_OUTPUT/PCoA"
#PCoA_INPUT="meta16S_OUTPUT/AlphaBetaDiversity/weighted_unifrac_otu_table.txt"
#mkdir $PCoA_DIR
#$qiime principal_coordinates.py -i $PCoA_INPUT -o $PCoA_DIR/beta_div_coords.txt
#$qiime make_2d_plots.py -i $PCoA_DIR/beta_div_coords.txt -m Map_file.txt -o $PCoA_DIR/

### Copy all meta16S_OUTPUT to ../meta16S_OUTPUT
sudo chown -R $username:$username *
rsync -avq meta16S_OUTPUT/* ../meta16S_OUTPUT

cd ..

echo "*** qiime work done! ***"

