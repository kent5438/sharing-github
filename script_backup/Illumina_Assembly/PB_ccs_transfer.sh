#! /bin/bash
### Create for NTUH BRCA Project
### Author: kentchen
### Date: 2018-07-13

ID=$1

usage() { echo "
Usage: $0 [-i <ID>]
-----
-i: The sequencing ID, (e.g. 001116)
-o: Output folder name
" && exit 1; }


while [ $# -gt 0 ]
do
OPT=$1

case $OPT in
    -h|--help)
    HELP=$1
    ((VAR++))
    ;;
    -i|--id)
    ID=$2
    shift
    ((VAR++))
    ;;
    -o|--output)
    OUT_DIR=$2
    shift
    ((VAR++))
    ;;
esac
shift
done

if [ ! ${ID} ]; then
	usage
elif [ ! ${OUT_DIR} ]; then
	usage
else
	echo "$(date +%T): % sh $0 -i $ID -o $OUT_DIR" >> CMD_history.txt
fi


JOBS_ROOT="/export/EC2480U/smrtlink_userdata_Ox/jobs_root"

if [ ! -e ${OUT_DIR} ]; then
	mkdir ${OUT_DIR}
fi
cd ${OUT_DIR}

PWD=`pwd`

ccs_fasta_dir="$JOBS_ROOT/001/${ID}/tasks/pbcoretools.tasks.bam2fasta_ccs-0/ccs.*.fasta"
ccs_fastq_dir="$JOBS_ROOT/001/${ID}/tasks/pbcoretools.tasks.bam2fastq_ccs-0/ccs.*.fastq"

rsync -av $ccs_fasta_dir $ccs_fastq_dir .

cd ..
