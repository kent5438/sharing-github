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

if [ ! -d "$MSFOLDER/data" ]; then
	mkdir -p "$MSFOLDER/data"
fi

cd "$MSFOLDER/data"

echo "*** Decompress..."

for i in `ls /mnt/NFS/EC2480U-P/scratch/QCResults/fastq/${MSID}/demux_cleanData/*R1*.fastq.gz | grep -v matched`; do zcat $i > $(basename $i | cut -d'_' -f1).R1.clean.fastq; done
for i in `ls /mnt/NFS/EC2480U-P/scratch/QCResults/fastq/${MSID}/demux_cleanData/*R2*.fastq.gz | grep -v matched`; do zcat $i > $(basename $i | cut -d'_' -f1).R2.clean.fastq; done

echo "*** Decompress complete *** "

perl /export/EC1680U/perl/bin/16S/CreateMappingFile.pl PE

echo "*** Successful <(_ _)> ***"

