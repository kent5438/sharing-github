#! /bin/bash
### Goals: (contigs.fasta) ---> (contigs.polished.fasta)

#now="[`date +"%F %T"`]"

prefix=$1
if [ -z ${prefix} ]; then
	echo "### ERROR: Give the prefix of sample name!!!
usage: $0 <prefix>";
	exit 1;
fi

check=`pwd | grep -c 'flye'`
if [ $check == 0 ]; then
	echo "### ERROR: Make sure you are located at '$(pwd)/flye'";
	exit 1;
fi

if [ ! -f "${prefix}.contigs.fasta" ]; then
	ln -s assembly.fasta ${prefix}.contigs.fasta
fi



### Main Pipeline

# LINKS
echo "### Run LINKS pipeline..." | ts '[%Y-%m-%d %H:%M:%S]'

if [ ! -f "../subreads/${prefix}.subreads.fastq" ]; then
	ls ../subreads/${prefix}.subreads.fasta > input.fofn
elif [ -f "../subreads/${prefix}.subreads.fastq" ]; then
	ls ../subreads/${prefix}.subreads.fastq > input.fofn
else
	echo "### ERROR: no subreads fastq/fasta found!"
fi

check_links=`ls | grep 'scaffolds.fa' | wc -l`

if [ ${check_links} == 0 ]; then
	/export/EC1680U/kentchen/opt/LINKS/releases/links_v1.8.7/LINKS -s input.fofn -f ${prefix}.contigs.fasta -b links
fi

# pbjelly
echo "### Run pbjelly pipeline..." | ts '[%Y-%m-%d %H:%M:%S]'

links_scaffolds=`ls *.scaffolds.fa`

if [ ! -f "scaffolds.fasta" ]; then
	ln -s ${links_scaffolds} scaffolds.fasta
fi

check_pbjelly=`ls | grep pbjelly_out | wc -l`
if [ ${check_pbjelly} == 0 ]; then
	sh /export/EC1680U/perl/bin/Pacbio/pbjelly.sh ${prefix}
fi


cd pbjelly_out

echo "Activate pbbioconda environment..." | ts '[%Y-%m-%d %H:%M:%S]'
. activate pbbioconda

echo "Generate fasta index..." | ts '[%Y-%m-%d %H:%M:%S]'
samtools faidx jelly.out.fasta

echo "Generate pbmm2 index..." | ts '[%Y-%m-%d %H:%M:%S]'
pbmm2 index jelly.out.fasta jelly.out.mmi

echo "Align raw subreads BAM to draft assembled contigs..." | ts '[%Y-%m-%d %H:%M:%S]'
bamCheck=`ls ../../rawData/*.bam | wc -l`
if [ $bamCheck != 1 ]; then
	echo "### ERROR: With / Without more than 1 raw subreads BAM in ../../rawData folder !!!"
	exit 1;
fi
bamFiles=`ls ../../rawData/*.bam`
pbmm2 align -j 32 --preset SUBREAD jelly.out.mmi ${bamFiles} | samtools sort -@ 32 -o ${prefix}.sorted.bam

echo "Generate BAM pbi index..." | ts '[%Y-%m-%d %H:%M:%S]'
pbindex ${prefix}.sorted.bam

echo "Polishing draft assembled contigs..." | ts '[%Y-%m-%d %H:%M:%S]'
arrow -j 32 ${prefix}.sorted.bam -r jelly.out.fasta -o ${prefix}.scaffolds.polished.fasta

cd ..

echo "Finish ~ <(_ _)>" | ts '[%Y-%m-%d %H:%M:%S]'
