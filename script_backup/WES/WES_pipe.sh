#! /bin/bash
### Create for NTUH WLHu WES Project
### Author: kentchen
### Date: 2017-12-25 ver 1

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


CLEAN_DIR=/export/EC2480U/WES/${ID}/cleanData
OUT_DIR=${ID}

if [ ! -e ${OUT_DIR} ]; then
	mkdir ${OUT_DIR}
fi
cd ${OUT_DIR}

PWD=`pwd`

for PREFIX in `ls ${CLEAN_DIR} | grep -v md5 | grep -v unpaired | awk -F'_R' '{print $1}' | sort -u`; 
do
	mkdir -p ${PREFIX}/reads
	ln -sf $(ls ${CLEAN_DIR}/${PREFIX}*R1*.fastq.gz | grep -v unpaired) ${PWD}/${PREFIX}/reads/${PREFIX}.R1.fastq.gz
	ln -sf $(ls ${CLEAN_DIR}/${PREFIX}*R2*.fastq.gz | grep -v unpaired) ${PWD}/${PREFIX}/reads/${PREFIX}.R2.fastq.gz
	rsync -avq /export/EC1680U/Pipeline/WES/GATK* ${PREFIX}/
	echo "${PREFIX}" > ${PREFIX}/list.txt
	cd ${PREFIX}/
	nohup time ./GATK_WES.ntuh.pl -s list.txt > nohup.out 2>nohup.err < /dev/null &
	echo $! > save_pid.txt
	cd ..
done





