#! /bin/bash

### only for 16S fastq so far ###

cycle=$1

if [ -z ${cycle} ]; then
	echo "### ERROR: give the cycle number!!!
usage: $0 <cycle>"
fi

. activate iDGE-20181208

if [ ! -e "fastq-stats-k" ]; then
	mkdir fastq-stats-k
fi

for i in `ls *R1*fastq*`; do 
	fastq-stats-k -D -c ${cycle} ${i} > fastq-stats-k/${i}-stats.txt; 
done

for i in `ls *R2*fastq*`; do
	fastq-stats-k -D -c ${cycle} ${i} > fastq-stats-k/${i}-stats.txt;
done

/mnt/NFS/EC2480U-P/scratch/iqc-pipeline/bin/prepareFastqStat.py fastq-stats-k

multiqc -c /mnt/NFS/EC2480U-P/scratch/iqc-pipeline/assets/multiqc_config.yaml -o multiqc_results fastq-stats-k
