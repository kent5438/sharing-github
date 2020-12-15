#! /bin/bash
### Create for NTUH BRCA Project
### Author: kentchen
### Date: 2018-07-13

ID=$1

usage() { echo "
Usage: $0 [-i <ID>] [-o <OUTPUT>]
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
	echo "[$(date)]: % sh $0 -i $ID -o $OUT_DIR" >> CMD_history.txt
fi


JOBS_ROOT="/export/EC2480U/smrtlink_userdata_Ox/jobs_root"

if [ ! -e ${OUT_DIR} ]; then
	mkdir ${OUT_DIR}
fi
cd ${OUT_DIR}

PWD=`pwd`

gff_dir="/export/EC2480U/smrtlink_userdata_Ox/jobs_root/001/${ID}/tasks/pbcoretools.tasks.gather_gff-1/file.gff"
vcf_dir="/export/EC2480U/smrtlink_userdata_Ox/jobs_root/001/${ID}/tasks/pbcoretools.tasks.gather_vcf-1/file.vcf"
bed_dir="/export/EC2480U/smrtlink_userdata_Ox/jobs_root/001/${ID}/tasks/genomic_consensus.tasks.gff2bed-0/output.bed"
aln_summary_gff="/export/EC2480U/smrtlink_userdata_Ox/jobs_root/001/${ID}/tasks/pbreports.tasks.summarize_coverage-0/alignment_summary.gff"

rsync -av $gff_dir $vcf_dir $bed_dir $aln_summary_gff .

grep "#" file.vcf > header

grep -v "#" file.vcf | grep BRCA1 | awk -F "\t" '{print "chr17\t"$2+43043552-1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9}' > BRCA1
grep -v "#" file.vcf | grep BRCA2 | awk -F "\t" '{print "chr13\t"$2+32314954-1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9}' > BRCA2

cat header BRCA2 BRCA1 > variant_calls.hg38.raw.vcf

rm -f header BRCA2 BRCA1

vep -i variant_calls.hg38.raw.vcf --cache --everything -o variant_calls.hg38.vep.vcf --dir /export/EC1680U/software/ensembl-vep/.vep --pick --pick_order appris --vcf --fork 30 --offline --species homo_sapiens --cache_version 92 --assembly GRCh38 --fasta /export/EC1680U/software/ensembl-vep/.vep/homo_sapiens/92_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

/export/EC1680U/perl/bin/WES/vcf2txt.v2.pl variant_calls.hg38.vep.vcf > variant_calls.hg38.vep.vcf.txt 2>/dev/null

/export/EC1680U/perl/bin/tab2xlsx.pl variant_calls.hg38.vep.vcf.txt variant_calls.hg38.vep.vcf.xlsx

rm -f variant_calls.hg38.vep.vcf.txt

cd ..
