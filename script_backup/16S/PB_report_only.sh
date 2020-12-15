#! /bin/bash
DIR=`pwd`
MSID=$1
MSFOLDER=$MSID.16S

sample_num=`cat SampleName_16S.txt | wc -l`

echo "[file]" > PB_report-16s.cfg
echo "BarcodeList=/export/EC1680U/kentchen/16S/barcode_list/${MSID}-16S-barcode_1.txt" >> PB_report-16s.cfg
echo "ReadStat=${DIR}/${sample_num}samples_readStat.txt" >> PB_report-16s.cfg
echo "OtuStat=${DIR}/meta16S_OUTPUT/OTUs_eachSample/${sample_num}samples.OTUstat.txt" >> PB_report-16s.cfg
echo "TaxonomyStat=${DIR}/meta16S_OUTPUT/OTUs_eachSample/${sample_num}samples.TAXONOMYstat.txt" >> PB_report-16s.cfg
echo "" >> PB_report-16s.cfg

echo "[template]" >> PB_report-16s.cfg
echo "Barcode=/export/EC1680U/kentchen/Template/PB_16S/barcodes_template" >> PB_report-16s.cfg
echo "Index=/export/EC1680U/kentchen/Template/PB_16S/index_template" >> PB_report-16s.cfg
echo "ReadStat=/export/EC1680U/kentchen/Template/PB_16S/readstat_template" >> PB_report-16s.cfg
echo "OtuStat=/export/EC1680U/kentchen/Template/PB_16S/otustat_template" >> PB_report-16s.cfg
echo "TaxonomyStat=/export/EC1680U/kentchen/Template/PB_16S/taxonomystat_template" >> PB_report-16s.cfg
echo "" >> PB_report-16s.cfg

echo "[parameter]" >> PB_report-16s.cfg
echo "SampleSize=$sample_num" >> PB_report-16s.cfg
echo "" >> PB_report-16s.cfg

ln -s ../QC.pl .
ln -s ../PB_report-16s.pl .

perl PB_report-16s.pl $MSID
perl QC.pl

mv -f report/PCA.pdf report/files
mv -f report/otu.table.forPCA.txt report/files

echo "### Complete!"

