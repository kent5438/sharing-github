#! /bin/bash

cycle=$1

if [ -z ${cycle} ]; then
	echo "### ERROR: give the cycle number!!!
usage: $0 <cycle>"
fi

. activate iDGE-20181208

cd reads

if [ ! -e "fastq-stats-k" ]; then
	mkdir fastq-stats-k
	
	for i in `ls *R1*fastq* | cut -d'_' -f1`; do
    	fastq-stats-k -D -c ${cycle} ${i}_R1.clean.fastq.gz > fastq-stats-k/${i}_R1.clean.fastq-stats.txt;
	done

	for i in `ls *R2*fastq* | cut -d'_' -f1`; do
    	fastq-stats-k -D -c ${cycle} ${i}_R2.clean.fastq.gz > fastq-stats-k/${i}_R2.clean.fastq-stats.txt;
	done
fi

/mnt/NFS/EC2480U-P/scratch/iqc-pipeline/bin/prepareFastqStat.py fastq-stats-k

cd ..

multiqc -c /export/EC1680U/perl/bin/Trinity/multiqc_config.yaml -o . reads/fastq-stats-k fastqc_result -f

