#! /bin/bash

DIR=`pwd`
echo "*** PWD: $DIR"

### Decompress CleanData & Run demultiplex program

MSID=$1
MSFOLDER=$MSID.ITS

if [ -z $MSID ]; then
	echo "### ERROR: Provide MiSeq project ID!"; 
	exit 1
fi

echo "*** CleanData path:  /export/EC2480U/Miseq/$MSID/CleanData"
rawdata_path="/export/EC2480U/Miseq/$MSID"

if [ ! -d $MSFOLDER ]; then
	mkdir $MSFOLDER
fi
cd $DIR/$MSID/CleanData

echo "*** Decompress..."
NUM_FQ_GZ=$(ls $DIR/$MSID/CleanData/*.fastq.gz | wc -l)
NUM_FQ=$(ls $DIR/$MSFOLDER/*.fastq 2>/dev/null | wc -l)
if [ "$NUM_FQ_GZ" -eq "$NUM_FQ" ]; then
	echo "### WARNING: $DIR/$MSID/CleanData/*.FASTQ.GZ have been DECOMPRESSED!"
else

echo "* fastq.gz unzip work..."
for i in $(ls *_001.clean.fastq.gz); 
do 
pigz -cd -p 32 ${MSID}-ITS_S*_L001_R1_001.clean.fastq.gz > "$DIR/$MSFOLDER/${MSID}-ITS.R1.clean.fastq"
pigz -cd -p 32 ${MSID}-ITS_S*_L001_R2_001.clean.fastq.gz > "$DIR/$MSFOLDER/${MSID}-ITS.R2.clean.fastq"
done | sort | uniq
echo "*** Decompress complete *** "

fi

cd $DIR/$MSFOLDER

if [ ! -e "$DIR/barcode_list/$MSID-ITS-barcode_1.txt" ]; then
touch "$DIR/barcode_list/$MSID-ITS-barcode_1.txt"
echo "### ERROR: Please complete barcode list file!";

exit 1

else
rsync -av "$DIR/barcode_list/$MSID-ITS-barcode_1.txt" "$DIR/$MSFOLDER"

fi

cd $DIR/$MSFOLDER
for i in $(ls *.fastq | cut -f 1 -d "." | sort | uniq) ; 
do echo "perl demultiplexFASTQ_byBarcodes.pl ${i}.R1.clean.fastq ${i}.R2.clean.fastq ${MSID}-ITS-barcode_1.txt _demultiplexed_${i} _baseSpace_${i} 2>&1 | tee -a demultiplexFASTQ-${i}.pl.log" ;


perl ../demultiplexFASTQ_byBarcodes.pl ${i}.R1.clean.fastq ${i}.R2.clean.fastq ${MSID}-ITS-barcode_1.txt _demultiplexed_${i} _baseSpace_${i} 2>&1 | tee -a demultiplexFASTQ-${i}.pl.log
done

### OTU process
cd $DIR/$MSFOLDER
mkdir 0_cleanFastq
ln -s /export/EC1680U/kentchen/ITS/$MSFOLDER/_demultiplexed_$MSID-ITS/*.fastq 0_cleanFastq 
cd 0_cleanFastq
rm -rf Undetermined*
ls *.R1.*fastq | cut -d "." -f 1 | sort -n > ../SampleName_ITS.txt
cd ..
ln -s $DIR/MiSeqITSanalysis_pipeline.v3v4.pl .
ln -s $DIR/V3FV4R_primer_V3V4.alignseq .

date

perl MiSeqITSanalysis_pipeline.v3v4.pl -p V3FV4R_primer_V3V4.alignseq -s SampleName_ITS.txt -r Y -o Y -i Y 2>&1 | tee -a MiSeqITSanalysis_pipeline.v3v4.MS16121.log

######################################################################
#  go to qiime VM and execute qiime_ITS_statistics_V3V4.kentchen.sh  #
######################################################################

echo "*** QIIME work in VM..."
ssh -t -t qiime@192.168.0.38<<EOF

cd /home/qiime/Desktop/OTU_analysis
mkdir $MSID; cd $MSID
scp kentchen@192.168.1.184:/export/EC1680U/kentchen/ITS/$MSFOLDER/6_biom/*.biom .
ln -s /home/qiime/Desktop/OTU_analysis/kentchen/qiime_ITS_statistics_V3V4.kentchen.sh .
./qiime_ITS_statistics_V3V4.kentchen.sh $MSID 2>&1 | tee qiime_process.log

exit

EOF

##########
# Report #
##########

cd $DIR/$MSFOLDER
echo "[file]" > report-ITS.cfg
echo "BarcodeList=/export/EC1680U/kentchen/ITS/$MSFOLDER/$MSID-ITS-barcode_1.txt" >> report-ITS.cfg
echo "ReadStat=/export/EC1680U/kentchen/ITS/$MSFOLDER/$(basename *samples_readStat.txt)" >> report-ITS.cfg
echo "OtuStat=/export/EC1680U/kentchen/ITS/$MSFOLDER/metaITS_OUTPUT/OTUs_eachSample/$(basename metaITS_OUTPUT/OTUs_eachSample/*samples.OTUstat.txt)" >> report-ITS.cfg
echo "TaxonomyStat=/export/EC1680U/kentchen/ITS/$MSFOLDER/metaITS_OUTPUT/OTUs_eachSample/$(basename metaITS_OUTPUT/OTUs_eachSample/*samples.TAXONOMYstat.txt)" >> report-ITS.cfg
echo "" >> report-ITS.cfg

echo "[template]" >> report-ITS.cfg
echo "Barcode=/export/EC1680U/kentchen/Template/ITS/barcodes_template" >> report-ITS.cfg
echo "Index=/export/EC1680U/kentchen/Template/ITS/index_template" >> report-ITS.cfg
echo "ReadStat=/export/EC1680U/kentchen/Template/ITS/readstat_template" >> report-ITS.cfg
echo "OtuStat=/export/EC1680U/kentchen/Template/ITS/otustat_template" >> report-ITS.cfg
echo "TaxonomyStat=/export/EC1680U/kentchen/Template/ITS/taxonomystat_template" >> report-ITS.cfg
echo "" >> report-ITS.cfg

echo "[parameter]" >> report-ITS.cfg
echo "SampleSize=$(cat SampleName_ITS.txt | wc -l)" >> report-ITS.cfg
echo "" >> report-ITS.cfg

if [ $(cat SampleName_ITS.txt | wc -l) -gt 2 ]; then
        Rscript /export/EC1680U/kentchen/ITS/PCA_plot.R $MSID
fi

ln -s $DIR/report-ITS.pl .
./report-ITS.pl $MSID

mv -f $DIR/$MSFOLDER/report/PCA.pdf $DIR/$MSFOLDER/report/files
rm -rf $DIR/$MSFOLDER/report/$MSID.otu.table.forPCA.txt

echo "### Please upload the $MSID.otu.table.forPCA.txt to the http://biit.cs.ut.ee/clustvis/ for PCA analysis."
echo "### zip -r $MSID report"
echo "### rsync -av $MSID.zip /export/EC1680U/USERDownLoad/<username>/ITS"


