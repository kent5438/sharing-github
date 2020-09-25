#! /bin/bash

DIR=`pwd`
echo "*** PWD: $DIR"

MSID=$1
MSFOLDER=$MSID.16S

if [ -z $MSID ]; then
	echo "### ERROR: Provide MiSeq project ID!"; 
	exit 1
fi


### OTU process
cd $DIR/$MSFOLDER
mkdir 0_cleanFastq
ln -s /export/EC1680U/kentchen/16S/$MSFOLDER/_demultiplexed_$MSID-16S/*.fastq 0_cleanFastq 
cd 0_cleanFastq
rm -rf Undetermined*
ls *.R1.*fastq | cut -d "." -f 1 | sort -n > ../SampleName_16S.txt
cd ..
ln -s $DIR/MiSeq16Sanalysis_pipeline.v3v4.pl .
ln -s $DIR/V3FV4R_primer_V3V4.alignseq .

date

perl MiSeq16Sanalysis_pipeline.v3v4.pl -p V3FV4R_primer_V3V4.alignseq -s SampleName_16S.txt -r Y -o Y -i Y 2>&1 | tee -a MiSeq16Sanalysis_pipeline.v3v4.MS16121.log

######################################################################
#  go to qiime VM and execute qiime_16S_statistics_V3V4.kentchen.sh  #
######################################################################

echo "*** QIIME work in VM..."
ssh -t -t qiime@192.168.0.38<<EOF

cd /home/qiime/Desktop/OTU_analysis
mkdir $MSID; cd $MSID
scp kentchen@192.168.1.184:/export/EC1680U/kentchen/16S/$MSFOLDER/6_biom/*.biom .
ln -s /home/qiime/Desktop/OTU_analysis/kentchen/qiime_16S_statistics_V3V4.kentchen.sh .
./qiime_16S_statistics_V3V4.kentchen.sh $MSID 2>&1 | tee qiime_process.log

exit

EOF

##########
# Report #
##########

cd $DIR/$MSFOLDER
echo "[file]" > report-16s.cfg
echo "BarcodeList=/export/EC1680U/kentchen/16S/$MSFOLDER/$MSID-16S-barcode_1.txt" >> report-16s.cfg
echo "ReadStat=/export/EC1680U/kentchen/16S/$MSFOLDER/$(basename *samples_readStat.txt)" >> report-16s.cfg
echo "OtuStat=/export/EC1680U/kentchen/16S/$MSFOLDER/meta16S_OUTPUT/OTUs_eachSample/$(basename meta16S_OUTPUT/OTUs_eachSample/*samples.OTUstat.txt)" >> report-16s.cfg
echo "TaxonomyStat=/export/EC1680U/kentchen/16S/$MSFOLDER/meta16S_OUTPUT/OTUs_eachSample/$(basename meta16S_OUTPUT/OTUs_eachSample/*samples.TAXONOMYstat.txt)" >> report-16s.cfg
echo "" >> report-16s.cfg

echo "[template]" >> report-16s.cfg
echo "Barcode=/export/EC1680U/kentchen/Template/16S/barcodes_template" >> report-16s.cfg
echo "Index=/export/EC1680U/kentchen/Template/16S/index_template" >> report-16s.cfg
echo "ReadStat=/export/EC1680U/kentchen/Template/16S/readstat_template" >> report-16s.cfg
echo "OtuStat=/export/EC1680U/kentchen/Template/16S/otustat_template" >> report-16s.cfg
echo "TaxonomyStat=/export/EC1680U/kentchen/Template/16S/taxonomystat_template" >> report-16s.cfg
echo "" >> report-16s.cfg

echo "[parameter]" >> report-16s.cfg
echo "SampleSize=$(cat SampleName_16S.txt | wc -l)" >> report-16s.cfg
echo "" >> report-16s.cfg

if [ $(cat SampleName_16S.txt | wc -l) -gt 2 ]; then
        Rscript /export/EC1680U/kentchen/16S/PCA_plot.R $MSID
fi

ln -s $DIR/report-16s.pl .
./report-16s.pl $MSID

mv -f $DIR/$MSFOLDER/report/PCA.pdf $DIR/$MSFOLDER/report/files
rm -rf $DIR/$MSFOLDER/report/$MSID.otu.table.forPCA.txt

echo "### Please upload the $MSID.otu.table.forPCA.txt to the http://biit.cs.ut.ee/clustvis/ for PCA analysis."
echo "### zip -r $MSID report"
echo "### rsync -av $MSID.zip /export/EC1680U/USERDownLoad/<username>/16S"
