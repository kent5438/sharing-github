#! /bin/bash
### Create for merge multiple fastq in different lane
### Author: kentchen
### Date: 2019-05-29

ID=$1

usage() { echo "
Usage: $0 [-i <ID>]
-----
-i: The sequencing ID, (e.g. NP17027)
" && exit 1; }


while getopts ":i:" o; do
	case "${o}" in
		i)
			ID=${OPTARG}
			;;
		*)
			usage
			;;
	esac
done
shift $((OPTIND-1))

if [ -z "${ID}" ]; then
	usage
else
	echo "$(date +%T): % sh $0 -i $ID" >> CMD_history.txt
fi


#CLEAN_DIR="/export/EC2480U/QCResults/fastq/${ID}/cleanData"
CLEAN_DIR="/mnt/NFS/EC2480U-P/scratch/QCResults/fastq/${ID}/cleanData"
OUT_DIR=${ID}

if [ ! -e "${OUT_DIR}/reads" ]; then
	mkdir -p ${OUT_DIR}/reads
fi
cd ${OUT_DIR}

PWD=`pwd`

echo "*** Start merge fastq..."

for PREFIX in `ls ${CLEAN_DIR}/*_L1_R1.clean.fastq.gz | awk -F'_' '{print $1}' | awk -F'cleanData/' '{print $2}'`; 
do
#	ln -sf $(ls ${CLEAN_DIR}/${PREFIX}*_R1_001*.fastq.gz | grep -v unpaired) ${PWD}/${PREFIX}/reads/${PREFIX}.R1.fastq.gz
	cat ${CLEAN_DIR}/${PREFIX}_*_L?_R1.clean.fastq.gz > reads/${PREFIX}_R1.clean.fastq.gz
#	ln -sf $(ls ${CLEAN_DIR}/${PREFIX}*_R2_001*.fastq.gz | grep -v unpaired) ${PWD}/${PREFIX}/reads/${PREFIX}.R2.fastq.gz
	cat ${CLEAN_DIR}/${PREFIX}_*_L?_R2.clean.fastq.gz > reads/${PREFIX}_R2.clean.fastq.gz
done
cd ..

echo "
*** Complete!
"



