#! /bin/bash

ID=`id -u -n`
PWD=`pwd`

OPT=$1 #{PB/Illu}

exist=`ls | grep 'D.1_taxa_summary' | wc -l`
thisFolder=`basename $PWD | grep 'OUTPUT' | wc -l`

if [ $exist -eq 0 ]; then
  echo 'Wrong folder!'
  exit
fi
if [ $thisFolder -eq 0 ]; then
  echo 'Wrong folder!'
  exit
fi

sudo chown -R ${ID}:users *

#if [ $exist -eq 0 ]; then
#	echo "### ERROR: this program should run on the project root folder!"
#	exit 1
#fi
#cd OUTPUT

if [ ! -d "report" ]; then
	mkdir report
fi

#if [ -d "E.3_LEfSe" ]; then
#	cd "E.3_LEfSe"
#	perl /export/EC1680U/perl/bin/16S/LEfSe_spec_charact.pl
#	cd ..
#fi

if [ "${OPT}" = "Illu" ]; then
	if [ ! -d B.5_OTU_cluster/97 ]; then
		ln -s $PWD/../OUTPUT/B.5_OTU_cluster/97 $PWD/B.5_OTU_cluster/97
		ln -s $PWD/../OUTPUT/B.1_QC $PWD/B.1_QC
		ln -s $PWD/../OUTPUT/B.2_join_fastq $PWD/B.2_join_fastq
		ln -s $PWD/../OUTPUT/B.3_demultiplex_quilty_filter $PWD/B.3_demultiplex_quilty_filter
		ln -s $PWD/../OUTPUT/B.4_removing_chimeras $PWD/B.4_removing_chimeras
	fi


	sudo docker run --rm  -v /etc/localtime:/etc/localtime:ro -w `pwd` -v /export/EC1680U/:/export/EC1680U/ -v /export/EC1680U/DataBase/16S:/references -v `pwd`:`pwd` kentchendocker/g16s_pipeline:latest /export/EC1680U/daniel/16s_template/plot_report.Illumina.rb 19 'report_html' '.' 'mothur_without_qiime' '/export/EC1680U/daniel/16s_template/template'

sudo chown -R $ID:users report_html

# Generate table from Clean -> NoPrimers -> Joined -> NoChimeras
	paste \
	<(cat B.1_QC/multiqc-rawdata/multiqc_data/multiqc_general_stats.txt | grep '\.R1' | sed 's/\.R1//' | cut -f1,5 | sort -k1,1 -V) \
	<(cat B.1_QC/multiqc-trimmed/multiqc_data/multiqc_general_stats.txt | grep '\.R1' | sed 's/\.R1//' | cut -f1,5 | sort -k1,1 -V | cut -f2) \
	<(for i in `ls B.2_join_fastq/*.fastq | sort -V`; do echo $(cat ${i} | wc -l)/4|bc; done) \
	<(for i in `ls B.4_removing_chimeras/*.trimmed.fa | sort -V`; do grep -c '>' ${i}; done) \
	<(for i in `ls ../OUTPUT/B_statistic_analysis/*_sample.txt | sort -V`; do cat ${i} | grep -v SampleID | cut -f2,3,4,5; done) \
	| sed 's/\.0//g' | perl -pe 's/^/SampleID\tClean\tNoPrimers\tJoined\tNoChimeras\tTotal_Tags\tUnique_Tags\tTaxon_Tags\tOTUs\n/ if $.==1' > B_statistic_analysis/reads-tags_summary.txt

        perl /export/EC1680U/perl/bin/16S/readTagSummary.pl

	rsync -avq B_statistic_analysis/reads-tags_summary.txt report_html/data/OTU_Statistic/

# Deriving Mapping.txt & BIOM for METAGENassist used
	mkdir report_html/data/METAGENassist
	output_underline_check=$(pwd | grep OUTPUT_Merge_ | wc -l)
	if [ $output_underline_check == 1 ]; then
		mapping_CompareId=$(basename `pwd` | sed 's/OUTPUT_Merge_//')
		mapping_txt="Mapping_${mapping_CompareId}.txt"
	else
		mapping_txt="Mapping.txt"
	fi

	# modified mapping.txt to mapping.csv (必須是csv結尾)
	cut -f1,4 ../data/`echo $mapping_txt` | sed 's/Description/Group/' | sed s/^\#// | sed 's/\t/,/' > report_html/data/METAGENassist/MetaData.csv
	rsync -avq B.5_OTU_cluster/otu_table.taxID.filter.dense.dedup.biom report_html/data/METAGENassist

###
        /usr/bin/perl /export/EC1680U/perl/bin/16S/16sFullTable.pl 
	
	rsync -avq E.1_PCA_PCoA/weighted_unifrac_pc_2D_PCoA_plots.png report_html/data/PCoA_weighted
	rsync -avq F.1_NMDS/PCoA_weighted_nmds.2d.* report_html/data/NMDS/Weighted
	rsync -avq /export/EC1680U/kentchen/16S/g16S/rmds/* report
	rsync -avq report_html/data  report
	# daniel add; kent modified if control
	if [ -d F.5_Tax4Fun ]; then
		rsync -avq F.5_Tax4Fun report/data
	fi
        sed 's/^/\"/' report/data/AlphaIndicesTable.txt | sed 's/$/\"/' | sed 's/\t/\",\"/g' > report/data/AlphaIndicesTable.csv
	cd report
	#sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` kendoit/g16s:0.4 \
	# daniel add
	sudo docker run --rm -v `pwd`:`pwd` -w `pwd` kendoit/g16s:0.5 \
		Rscript -e "rmarkdown::render_site()"
	cd ..
	sudo chown -R ${ID}:users report
	rsync -avq B.5_OTU_cluster/otu_table.taxID.filter.dense.dedup.biom report/_reports/data/
	rsync -avq B.5_OTU_cluster/otu_table.taxID.filter.dense.dedup.txt report/_reports/data/

### copy beta diversity result
	if [ ! -d report/_reports/data/BetaDiversity ]; then
		mkdir report/_reports/data/BetaDiversity
		rsync -avq E.1_PCA_PCoA/*_dm.txt report/_reports/data/BetaDiversity
	fi
###

### parse taxonomy info of rep fasta for user (powered by JSYF)
	perl /export/EC1680U/perl/bin/16S/otuCluster.pl
	rm -f report/_reports/advanced.Rmd.20190429
###

elif [ "${OPT}" = "PB" ]; then
        if [ ! -d B.5_OTU_cluster/97 ]; then
                ln -s $PWD/../OUTPUT/B.5_OTU_cluster/97 $PWD/B.5_OTU_cluster/97

                ln -s $PWD/../OUTPUT/B.1_QC $PWD/B.1_QC
#                ln -s $PWD/../OUTPUT/B.2_join_fastq $PWD/B.2_join_fastq
                ln -s $PWD/../OUTPUT/B.3_demultiplex_quilty_filter $PWD/B.3_demultiplex_quilty_filter
                ln -s $PWD/../OUTPUT/B.4_removing_chimeras $PWD/B.4_removing_chimeras
        fi

	sudo docker run -it --rm  -v /etc/localtime:/etc/localtime:ro -w `pwd` -v /export/EC1680U/:/export/EC1680U/ -v /export/EC1680U/DataBase/16S:/references -v `pwd`:`pwd` kentchendocker/g16s_pipeline:latest /export/EC1680U/daniel/16s_template/plot_report.PacBio.rb 19 'report_html' '.' 'mothur_without_qiime' '/export/EC1680U/daniel/16s_template/template'

    sudo chown -R $ID:users report_html

# Generate table from Clean -> NoPrimers -> Joined -> NoChimeras
	paste \
	<(cat B.1_QC/multiqc-rawdata/multiqc_data/multiqc_general_stats.txt | grep -v '^Sample' | cut -f1,5 | sort -k1,1 -V) \
	<(cat B.1_QC/multiqc-trimmed/multiqc_data/multiqc_general_stats.txt | grep -v '^Sample' | cut -f1,5 | sort -k1,1 -V | cut -f2) \
	<(for i in `ls B.4_removing_chimeras/*.trimmed.fa | sort -V`; do grep -c '>' ${i}; done) \
	<(for i in `ls B_statistic_analysis/*_sample.txt | sort -V`; do cat ${i} | grep -v SampleID | cut -f2,3,4,5; done) \
	| sed 's/\.0//g' | perl -pe 's/^/SampleID\tClean\tNoPrimers\tNoChimeras\tTotal_Tags\tUnique_Tags\tTaxon_Tags\tOTUs\n/ if $.==1' > B_statistic_analysis/reads-tags_summary.txt

	perl /export/EC1680U/perl/bin/16S/readTagSummary.PB.pl

	rsync -avq B_statistic_analysis/reads-tags_summary.txt report_html/data/OTU_Statistic/

# Deriving MetaData.csv & BIOM for METAGENassist used
    mkdir report_html/data/METAGENassist
    output_underline_check=$(pwd | grep OUTPUT_Merge_ | wc -l)
    if [ $output_underline_check == 1 ]; then
        mapping_CompareId=$(basename `pwd` | sed 's/OUTPUT_Merge_//')
        mapping_txt="Mapping_${mapping_CompareId}.txt"
    else
        mapping_txt="Mapping.txt"
    fi

# modified mapping.txt to mapping.csv (必須是csv結尾)
    cut -f1,4 ../data/`echo $mapping_txt` | sed 's/Description/Group/' | sed s/^\#// | sed 's/\t/,/' > report_html/data/METAGENassist/MetaData.csv
    rsync -avq B.5_OTU_cluster/otu_table.taxID.filter.dense.dedup.biom report_html/data/METAGENassist

###

    rsync -avq E.1_PCA_PCoA/weighted_unifrac_pc_2D_PCoA_plots.png report_html/data/PCoA_weighted
    rsync -avq F.1_NMDS/PCoA_weighted_nmds.2d.* report_html/data/NMDS/Weighted
    rsync -avq /export/EC1680U/kentchen/16S/g16S/rmds/* report
    rsync -avq report_html/data  report
    # daniel add; kent modified if control
    if [ -d F.5_Tax4Fun ]; then
        rsync -avq F.5_Tax4Fun report/data
    fi


	rsync -avq /export/EC1680U/kentchen/16S/g16S/rmds_pacbio/* report
	rsync -avq report_html/data report
	if [ ! -d B.5_OTU_cluster/97 ]; then
		ln -s $PWD/../OUTPUT/B.5_OTU_cluster/97 $PWD/B.5_OTU_cluster/97
	fi
        sed 's/^/\"/' report/data/AlphaIndicesTable.txt | sed 's/$/\"/' | sed 's/\t/\",\"/g' > report/data/AlphaIndicesTable.csv
	cd report
	sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` kendoit/g16s:0.5 \
		Rscript -e "rmarkdown::render_site()"
	cd ..
	sudo chown -R ${ID}:users report
	rsync -avq B.5_OTU_cluster/otu_table.taxID.filter.dense.dedup.biom report/_reports/data
	rsync -avq B.5_OTU_cluster/otu_table.taxID.filter.dense.dedup.txt report/_reports/data
### parse more detail rep fasta for user (powered by JSYF)
### copy beta diversity result
        if [ ! -d report/_reports/data/BetaDiversity ]; then
                mkdir report/_reports/data/BetaDiversity
                rsync -avq E.1_PCA_PCoA/*_dm.txt report/_reports/data/BetaDiversity
        fi
###
### parse taxonomy info of rep fasta for user (powered by JSYF)
    perl /export/EC1680U/perl/bin/16S/otuCluster.pl
    rm -f report/_reports/advanced.Rmd.20190429
###

else
	echo "### ERROR: please make sure the 1st ARGS should be [PB/Illu]"
	exit 1
fi
rsync -avq report_html/report.pdf report/_reports
rsync -avq /export/EC1680U/kentchen/Template/16S/METAGENassist_manual.pdf report/_reports/data/METAGENassist/
cd ..

