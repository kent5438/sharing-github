#! /bin/bash
DIR=`pwd`
MSID=$1
MSFOLDER=$MSID.ITS

sample_num=`cat SampleName_ITS.txt | wc -l`

echo "[file]" > report-ITS.cfg
echo "BarcodeList=/export/EC1680U/kentchen/ITS/barcode_list/$MSID-ITS-barcode_1.txt" >> report-ITS.cfg
echo "ReadStat=${DIR}/${sample_num}samples_readStat.txt" >> report-ITS.cfg
echo "OtuStat=${DIR}/metaITS_OUTPUT/OTUs_eachSample/${sample_num}samples.OTUstat.txt" >> report-ITS.cfg
echo "TaxonomyStat=${DIR}/metaITS_OUTPUT/OTUs_eachSample/${sample_num}samples.TAXONOMYstat.txt" >> report-ITS.cfg
echo "" >> report-ITS.cfg

echo "[template]" >> report-ITS.cfg
echo "Barcode=/export/EC1680U/kentchen/Template/ITS/barcodes_template" >> report-ITS.cfg
echo "Index=/export/EC1680U/kentchen/Template/ITS/index_template" >> report-ITS.cfg
echo "ReadStat=/export/EC1680U/kentchen/Template/ITS/readstat_template" >> report-ITS.cfg
echo "OtuStat=/export/EC1680U/kentchen/Template/ITS/otustat_template" >> report-ITS.cfg
echo "TaxonomyStat=/export/EC1680U/kentchen/Template/ITS/taxonomystat_template" >> report-ITS.cfg
echo "" >> report-ITS.cfg

echo "[parameter]" >> report-ITS.cfg
echo "SampleSize=$sample_num" >> report-ITS.cfg
echo "" >> report-ITS.cfg

ln -s /export/EC1680U/kentchen/ITS/QC.pl .
ln -s /export/EC1680U/kentchen/ITS/report-ITS.pl .

perl report-ITS.pl $MSID
perl QC.pl

mv -f report/PCA.pdf report/files
mv -f report/otu.table.forPCA.txt report/files

echo "### Complete!"
