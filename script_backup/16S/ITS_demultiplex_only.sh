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
if [ ! -e "$DIR/$MSFOLDER/${MSID}-ITS.*.clean.fastq" ]; then
for i in $(ls *_001.clean.fastq.gz); 
do 
pigz -cd -p 32 ${MSID}-ITS_S*_R1_001.clean.fastq.gz > "$DIR/$MSFOLDER/${MSID}-ITS.R1.clean.fastq"
pigz -cd -p 32 ${MSID}-ITS_S*_R2_001.clean.fastq.gz > "$DIR/$MSFOLDER/${MSID}-ITS.R2.clean.fastq"
done | sort | uniq
fi

echo "*** Decompress complete *** "

fi

cd $DIR/$MSFOLDER

if [ ! -e "$DIR/barcode_list/$MSID-ITS-barcode_1.txt" ]; then
touch "$DIR/barcode_list/$MSID-ITS-barcode_1.txt"
echo "### ERROR: Please complete barcode list file!";

#dualbarcodes	barcode1        barcode2	sampleName
dualbarcodes	GCGTAGTA	CTCTCTAT	1
dualbarcodes	CGGAGCCT	CTCTCTAT	2
dualbarcodes	TACGCTGC	CTCTCTAT	3
dualbarcodes	ATGCGCAG	CTCTCTAT	4
dualbarcodes	TAGCGCTC	CTCTCTAT	5
dualbarcodes	ACTGAGCG	CTCTCTAT	6

exit 1

else
rsync -av "$DIR/barcode_list/$MSID-ITS-barcode_1.txt" "$DIR/$MSFOLDER"

fi

cd $DIR/$MSFOLDER
for i in $(ls *.fastq | cut -f 1 -d "." | sort | uniq) ; 
do echo "perl demultiplexFASTQ_byBarcodes.pl ${i}.R1.clean.fastq ${i}.R2.clean.fastq ${MSID}-ITS-barcode_1.txt _demultiplexed_${i} _baseSpace_${i} 2>&1 | tee -a demultiplexFASTQ-${i}.pl.log" ;

perl ../demultiplexFASTQ_byBarcodes.pl ${i}.R1.clean.fastq ${i}.R2.clean.fastq ${MSID}-ITS-barcode_1.txt _demultiplexed_${i} _baseSpace_${i} 2>&1 | tee -a demultiplexFASTQ-${i}.pl.log
done

