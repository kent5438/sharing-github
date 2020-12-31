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

/opt/smrtlink/smrtcmds/bin/ccs --minPasses=1 --reportFile=ccs.report /export/EC2480U/smrtlink_userdata_Ox/jobs_root/001/${ID}/merged.dataset.xml ccs.bam
/opt/smrtlink/smrtcmds/bin/bam2fastq -u -o ccs ccs.bam
rm -f ccs.bam*

ngmlr -t 32 -r /export/EC1680U/kentchen/Customized/NTUH_BRCA-simon/Ref/BRCA1-2_NTUH-simon.fasta -q ccs.fastq  | samtools sort -@ 32 -o ${OUT_DIR}.sorted.bam -

sniffles -m ${OUT_DIR}.sorted.bam -b ${OUT_DIR}.bedpe -t 32; /export/EC1680U/perl/bin/tab2xlsx.pl ${OUT_DIR}.bedpe ${OUT_DIR}.bedpe.xlsx; sniffles -m ${OUT_DIR}.sorted.bam -v ${OUT_DIR}.vcf; /export/EC1680U/perl/bin/tab2xlsx.pl ${OUT_DIR}.vcf ${OUT_DIR}.vcf.xlsx

cd ..
