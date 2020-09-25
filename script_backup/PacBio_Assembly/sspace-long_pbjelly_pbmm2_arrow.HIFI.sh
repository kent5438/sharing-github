#! /bin/bash
### Goals: (assembly.fasta) ---> (contigs.polished.fasta)

#now="[`date +"%F %T"`]"

prefix=$1
if [ -z ${prefix} ]; then
	echo "### ERROR: Give the prefix of sample name!!!
usage: $0 <prefix>";
	exit 1;
fi

check=`pwd | grep -c 'hifi'`
if [ $check == 0 ]; then
	echo "### ERROR: Make sure you are located at '$(pwd)/hifi'";
	exit 1;
fi

if [ ! -f "${prefix}.contigs.fasta" ]; then
	ln -s hifiasm.p_ctg.fasta ${prefix}.contigs.fasta
fi



### Main Pipeline

echo "### Check scaffolds.polished.fasta existed or not..." | ts '[%Y-%m-%d %H:%M:%S]'
if [ -f "PacBio_scaffolder_results/pbjelly_out/${prefix}.scaffolds.polished.fasta" ]; then
	echo "### Hey! You already have finished this pipeline!"
	exit 1;
fi

### SSPACE-Long

echo "### Running SSPACE-Long scaffolding..." | ts '[%Y-%m-%d %H:%M:%S]'
subreads=`ls ../subreads/*.fastq`
if [ `ls ../subreads/*.fastq | wc -l` != 1 ]; then
    echo "### ERROR: Without / With more than 1 subreads fastq in ../subreads folder !!!"
    exit 1;
fi

if [ ! -f "PacBio_scaffolder_results/scaffolds.fasta" ]; then
	perl /export/EC1680U/kentchen/opt/SSPACE-LongRead_v1-1/SSPACE-LongRead_minimap2.pl -c ${prefix}.contigs.fasta -p ${subreads} -t 32
fi

cd PacBio_scaffolder_results


### pbjelly

echo "### Run pbjelly pipeline..." | ts '[%Y-%m-%d %H:%M:%S]'

if [ ! -f "scaffolds.fasta" ]; then
    echo "### ERROR: make sure you have complete SSPACE-Long pipeline"
	exit 1;
fi

check_pbjelly=`ls | grep pbjelly_out | wc -l`
if [ ${check_pbjelly} == 0 ]; then
    sh /export/EC1680U/perl/bin/Pacbio/pbjelly.sh ${prefix}
fi


cd pbjelly_out


### Polishing

echo "### Activate pbbioconda environment..." | ts '[%Y-%m-%d %H:%M:%S]'
. activate pbbioconda

echo "### Generate fasta index..." | ts '[%Y-%m-%d %H:%M:%S]'
samtools faidx jelly.out.fasta

echo "### Generate pbmm2 index..." | ts '[%Y-%m-%d %H:%M:%S]'
pbmm2 index jelly.out.fasta jelly.out.mmi

echo "### Align raw subreads BAM to draft assembled contigs..." | ts '[%Y-%m-%d %H:%M:%S]'
bamCheck=`ls ../../../rawData/*.bam | wc -l`
if [ $bamCheck != 1 ]; then
	echo "### ERROR: Without / With more than 1 raw subreads BAM in ../rawData folder !!!"
	exit 1;
fi
bamFiles=`ls ../../../rawData/*.bam`
pbmm2 align --sort -j 32 --preset SUBREAD jelly.out.mmi ${bamFiles} ${prefix}.sorted.bam

echo "### Generate BAM pbi index..." | ts '[%Y-%m-%d %H:%M:%S]'
pbindex ${prefix}.sorted.bam

echo "### Polishing draft assembled contigs..." | ts '[%Y-%m-%d %H:%M:%S]'
arrow -j 32 ${prefix}.sorted.bam -r jelly.out.fasta -o ${prefix}.scaffolds.polished.fasta

echo "### pbmm2 + arrow Finish ~ <(_ _)>" | ts '[%Y-%m-%d %H:%M:%S]'

cd ..
