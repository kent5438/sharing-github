#! /bin/bash
DIR=`pwd`
echo "*** PWD: $DIR"

### Decompress CleanData & Run demultiplex program

MSID=$1
#MSID_nodash=`echo ${MSID//-/ } | cut -d' ' -f1`


MSFOLDER=$DIR/${MSID}.16S
#MiseqClean="/export/EC2480U/QCResults/fastq/$MSID/cleanData"
MiseqClean="/mnt/NFS/EC2480U-P/scratch/QCResults/fastq/$MSID/cleanData"

if [ -z $MSID ]; then
	echo "### ERROR: Provide MiSeq project ID!"; 
	exit 1
fi

if [ ! -d $MSFOLDER ]; then
	mkdir $MSFOLDER
fi


echo "*** Decompress..."

if [ ! -e "$MSFOLDER/${MSID}.*.clean.fastq" ]; then
	/export/EC1680U/software/anaconda2/bin/pigz -cd -p 32 `ls ${MiseqClean}/${MSID}-pool*R1*.clean.fastq.gz | grep -v unpaired` > "$MSFOLDER/${MSID}.R1.clean.fastq"
	/export/EC1680U/software/anaconda2/bin/pigz -cd -p 32 `ls ${MiseqClean}/${MSID}-pool*R2*.clean.fastq.gz | grep -v unpaired` > "$MSFOLDER/${MSID}.R2.clean.fastq"
fi;


echo "*** Decompress complete *** "

cd $MSFOLDER

if [ ! -e "$DIR/barcode_list/$MSID-16S-barcode_1.txt" ]; then
touch "$DIR/barcode_list/$MSID-16S-barcode_1.txt"
echo "### ERROR: Please complete barcode list file!";
exit 1

else
rsync -av "$DIR/barcode_list/$MSID-16S-barcode_1.txt" "$MSFOLDER"

fi

for i in $(ls ${MSID}.R*.clean.fastq | cut -f 1 -d "." | sort | uniq) ; 
do echo "perl /export/EC1680U/perl/bin/16S/demultiplexFASTQ_byBarcodes.pl ${i}.R1.clean.fastq ${i}.R2.clean.fastq ${MSID}-16S-barcode_1.txt _demultiplexed_${i} _baseSpace_${i} 2>&1 | tee -a demultiplexFASTQ-${i}.pl.log" ;

/export/EC1680U/kentchen/miniconda3/envs/genomics-16S/bin/perl /export/EC1680U/perl/bin/16S/demultiplexFASTQ_byBarcodes.pl ${i}.R1.clean.fastq ${i}.R2.clean.fastq ${MSID}-16S-barcode_1.txt _demultiplexed_${i} _baseSpace_${i} 2>&1 | tee -a demultiplexFASTQ-${i}.pl.log
done


mkdir 0_cleanFastq 

cd 0_cleanFastq
ln -s $(ls ../_demultiplexed_${MSID}/*.clean.fastq | grep -v Undetermined) .
cd ..

ln -s /export/EC1680U/perl/bin/16S/MiSeq16Sanalysis_pipeline.v3v4cutadapt.SILVA.pl .
ln -s ../V3FV4R_primer_V3V4.alignseq .
ls 0_cleanFastq/*.clean.fastq | cut -d'.' -f1 | cut -d'/' -f2 | sort -u > SampleName_16S.txt

echo ""
echo "perl MiSeq16Sanalysis_pipeline.v3v4cutadapt.SILVA.pl -s SampleName_16S.txt -p V3FV4R_primer_V3V4.alignseq -r Y -o Y -i N 2>&1 | tee mothur.log"
echo ""


