#! /bin/bash
DIR=`pwd`
MSID=$1
MSFOLDER=$MSID.16S

sample_num=`cat SampleName_16S.txt | wc -l`

echo "[file]" > report-16s.cfg
echo "BarcodeList=/export/EC1680U/kentchen/16S/barcode_list/$MSID-16S-barcode_1.txt" >> report-16s.cfg
echo "ReadStat=${DIR}/${sample_num}samples_readStat.txt" >> report-16s.cfg
echo "OtuStat=${DIR}/meta16S_OUTPUT/OTUs_eachSample/${sample_num}samples.OTUstat.txt" >> report-16s.cfg
echo "TaxonomyStat=${DIR}/meta16S_OUTPUT/OTUs_eachSample/${sample_num}samples.TAXONOMYstat.txt" >> report-16s.cfg
echo "" >> report-16s.cfg

echo "[template]" >> report-16s.cfg
echo "Barcode=/export/EC1680U/kentchen/Template/16S/barcodes_template" >> report-16s.cfg
echo "Index=/export/EC1680U/kentchen/Template/16S/index_template" >> report-16s.cfg
echo "ReadStat=/export/EC1680U/kentchen/Template/16S/readstat_template" >> report-16s.cfg
echo "OtuStat=/export/EC1680U/kentchen/Template/16S/otustat_template" >> report-16s.cfg
echo "TaxonomyStat=/export/EC1680U/kentchen/Template/16S/taxonomystat_template" >> report-16s.cfg
echo "" >> report-16s.cfg

echo "[parameter]" >> report-16s.cfg
echo "SampleSize=$sample_num" >> report-16s.cfg
echo "" >> report-16s.cfg

perl /export/EC1680U/kentchen/16S/report-16s.pl $MSID
perl /export/EC1680U/kentchen/16S/QC.pl

mv -f report/PCA.pdf report/files
mv -f report/otu_table.forPCA.txt report/files

echo "### Complete!"
