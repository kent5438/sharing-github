#! /bin/bash

DIR=`pwd`
MSID=$1
MSFOLDER=$1.16S

echo "*** QIIME work in VM..."
ssh -t -t qiime@192.168.0.38<<EOF

cd /home/qiime/Desktop/OTU_analysis
mkdir $MSID; cd $MSID
scp kentchen@192.168.1.184:/export/EC1680U/kentchen/16S/$MSFOLDER/6_biom/*.biom .
ln -s /home/qiime/Desktop/OTU_analysis/kentchen/qiime_16S_statistics_V3V4.kentchen.sh .
./qiime_16S_statistics_V3V4.kentchen.sh $MSID 2>&1 | tee qiime_process.log

exit
EOF

cd $DIR/$MSFOLDER
sample_num=`cat SampleName_16S.txt | wc -l`

echo "[file]" > report-16s.cfg
echo "BarcodeList=/export/EC1680U/kentchen/16S/$MSFOLDER/$MSID-16S-barcode_1.txt" >> report-16s.cfg
echo "ReadStat=/export/EC1680U/kentchen/16S/$MSFOLDER/${sample_num}samples_readStat.txt" >> report-16s.cfg
echo "OtuStat=/export/EC1680U/kentchen/16S/$MSFOLDER/meta16S_OUTPUT/OTUs_eachSample/meta16S_OUTPUT/OTUs_eachSample/${sample_num}samples.OTUstat.txt" >> report-16s.cfg
echo "TaxonomyStat=/export/EC1680U/kentchen/16S/$MSFOLDER/meta16S_OUTPUT/OTUs_eachSample/meta16S_OUTPUT/OTUs_eachSample/${sample_num}samples.TAXONOMYstat.txt" >> report-16s.cfg
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

ln -s $DIR/report-16s.pl .

if [ $(cat SampleName_16S.txt | wc -l) -gt 2 ]; then
        Rscript /export/EC1680U/kentchen/16S/PCA_plot.R $MSID
fi

echo "*** Generate report folder..."
./report-16s.pl

mv -f $DIR/$MSFOLDER/report/PCA.pdf $DIR/$MSFOLDER/report/files
rm -rf $DIR/$MSFOLDER/report/$MSID.otu.table.forPCA.txt

