#! /bin/bash
### Goals: (contigs.fasta) ---> (contigs.polished.fasta)

#now="[`date +"%F %T"`]"

prefix=$1
if [ -z ${prefix} ]; then
	echo "### ERROR: Give the prefix of sample name!!!
usage: $0 <prefix>";
	exit 1;
fi

check=`pwd | grep -c 'wtdbg2'`
if [ $check == 0 ]; then
	echo "### ERROR: Make sure you are located at '$(pwd)/wtdbg2'";
	exit 1;
fi

if [ ! -f "${prefix}.contigs.fasta" ]; then
	ln -s dbg.cns.ctg.fa ${prefix}.contigs.fasta
fi



# Main Pipeline

echo "Check contigs.polished.fasta existed or not..." | ts '[%Y-%m-%d %H:%M:%S]'
if [ -f "${prefix}.contigs.polished.fasta" ]; then
	echo "### Hey! You already have the ${prefix}.contigs.polished.fasta file!"
	exit 1;
fi

echo "Activate wtdbg2 environment..." | ts '[%Y-%m-%d %H:%M:%S]'
. activate pbbioconda

echo "Generate fasta index..." | ts '[%Y-%m-%d %H:%M:%S]'
samtools faidx ${prefix}.contigs.fasta

echo "Generate pbmm2 index..." | ts '[%Y-%m-%d %H:%M:%S]'
pbmm2 index ${prefix}.contigs.fasta ${prefix}.contigs.mmi

echo "Align raw subreads BAM to draft assembled contigs..." | ts '[%Y-%m-%d %H:%M:%S]'
bamCheck=`ls ../rawData/*.bam | wc -l`
if [ $bamCheck != 1 ]; then
	echo "### ERROR: Without / With more than 1 raw subreads BAM in ../rawData folder !!!"
	exit 1;
fi
bamFiles=`ls ../rawData/*.bam`
pbmm2 align --sort -j 32 --preset SUBREAD ${prefix}.contigs.mmi ${bamFiles} ${prefix}.sorted.bam

echo "Generate BAM pbi index..." | ts '[%Y-%m-%d %H:%M:%S]'
pbindex ${prefix}.sorted.bam

echo "Polishing draft assembled contigs..." | ts '[%Y-%m-%d %H:%M:%S]'
arrow -j 32 ${prefix}.sorted.bam -r ${prefix}.contigs.fasta -o ${prefix}.contigs.polished.fasta

echo "Finish ~ <(_ _)>" | ts '[%Y-%m-%d %H:%M:%S]'
