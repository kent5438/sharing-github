#! /bin/bash
DIR=`pwd`
echo "*** PWD: $DIR"

### Decompress CleanData & Run demultiplex program

MSID=$1
MSID_nodash=(${MSID//-/ })


MSFOLDER=$DIR/${MSID}.16S
MiseqClean="/export/EC2480U/QCResults/fastq/$MSID_nodash/cleanData"

#if [ -z $MSID ]; then
#	echo "### ERROR: Provide MiSeq project ID!"; 
#	exit 1
#fi

#if [ ! -d $MSFOLDER ]; then
#	mkdir $MSFOLDER
#fi


#echo "*** Decompress..."
#NUM_FQ_GZ=$(ls $MiseqClean/*.clean.fastq.gz | grep -v unpaired | wc -l)
#NUM_FQ=$(ls $MSFOLDER/*.fastq 2>/dev/null | wc -l)
#if [ "$NUM_FQ_GZ" -eq "$NUM_FQ" ]; then
#	echo "### WARNING: $DIR/$MSID/cleanData/*.FASTQ.GZ have been DECOMPRESSED!"
#else
#
#echo "* fastq.gz unzip work..."
#for i in $(ls ${MiseqClean}/*.clean.fastq.gz | grep -v unpaired | awk -F'_R' '{print $1}');
#do
#	if [ ! -e "$MSFOLDER/${MSID}.*.clean.fastq" ]; then
#		/export/EC1680U/software/anaconda2/bin/pigz -cd -p 32 `ls ${MiseqClean}/*_R1_00?.clean.fastq.gz | grep -v unpaired` > "$MSFOLDER/${MSID}.R1.clean.fastq"
#		/export/EC1680U/software/anaconda2/bin/pigz -cd -p 32 `ls ${MiseqClean}/*_R2_00?.clean.fastq.gz | grep -v unpaired` > "$MSFOLDER/${MSID}.R2.clean.fastq"
#	fi;
#done | sort | uniq


#echo "*** Decompress complete *** "

#fi
cd $MSFOLDER

if [ ! -e "$DIR/barcode_list/$MSID-16S-barcode_1.txt" ]; then
touch "$DIR/barcode_list/$MSID-16S-barcode_1.txt"
echo "### ERROR: Please complete barcode list file!";
exit 1

else
rsync -av "$DIR/barcode_list/$MSID-16S-barcode_1.txt" "$MSFOLDER"

fi

#for i in $(ls ${MSID}.R*.clean.fastq | cut -f 1 -d "." | sort | uniq) ; 
#do echo "perl demultiplexFASTQ_byBarcodes.pl ${i}.R1.clean.fastq ${i}.R2.clean.fastq ${MSID}-16S-barcode_1.txt _demultiplexed_${i} _baseSpace_${i} 2>&1 | tee -a demultiplexFASTQ-${i}.pl.log" ;

#perl ../demultiplexFASTQ_byBarcodes.pl ${i}.R1.clean.fastq ${i}.R2.clean.fastq ${MSID}-16S-barcode_1.txt _demultiplexed_${i} _baseSpace_${i} 2>&1 | tee -a demultiplexFASTQ-${i}.pl.log
#done


#mkdir 0_cleanFastq 

cd 0_cleanFastq
ln -s $(ls $MSFOLDER/_demultiplexed_${MSID}/*.clean.fastq | grep -v Undetermined) .
cd ..

ln -s /export/EC1680U/perl/bin/16S/MiSeq16Sanalysis_pipeline.v3v4cutadapt.SILVA.pl .
ln -s ../V3FV4R_primer_V3V4.alignseq .
ls 0_cleanFastq/*.clean.fastq | cut -d'.' -f1 | cut -d'/' -f2 | sort -u > SampleName_16S.txt

echo ""
echo "perl MiSeq16Sanalysis_pipeline.v3v4cutadapt.SILVA.pl -s SampleName_16S.txt -p V3FV4R_primer_V3V4.alignseq -r Y -o Y -i N 2>&1 | tee mothur.log"
echo ""


