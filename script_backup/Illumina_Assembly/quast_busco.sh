#! /bin/bash
### Goals: generate assembly summary for (contigs.polished.fasta)

### USAGE:
# % ./quast_busco.sh ${prefix} ${busco_db}

prefix=$1
if [ -z ${prefix} ]; then
    echo "### ERROR: Give the prefix of sample name!!!
usage: $0 <prefix> <busco_db>";
    exit 1;
fi

busco_db=$2
if [ -z ${busco_db} ]; then
    echo "### ERROR: Give the busco database full path!!!
usage: $0 <prefix> <busco_db>";
    exit 1;
fi

check=`pwd | grep -c 'flye'`
if [ $check == 0 ]; then
    echo "### ERROR: Make sure you are located at '$(pwd)/flye'";
    exit 1;
fi

if [ ! -f "${prefix}.contigs.polished.fasta" ]; then
	echo "### ERROR: please make sure you have completed 'pbmm2_arrow.flye.sh'"
	exit 1;
fi


# Main Pipeline
echo "### Assembly evaluation by QUAST..." | ts '[%Y-%m-%d %H:%M:%S]'
. activate quast
quast -t 32 ${prefix}.contigs.polished.fasta

echo "### Assembly evaluation by BUSCO..." | ts '[%Y-%m-%d %H:%M:%S]'
. activate busco
run_BUSCO.py -i ${i}.contigs.polished.fasta -o busco -m geno -l ${busco_db} -t /home/kentchen/tmp -f -c 32 --blast_single_core

echo "### Generate report folder..." | ts '[%Y-%m-%d %H:%M:%S]'
mkdir report
cd report
mkdir 01_QUAST 02_Assembly 03_BUSCO
rsync -avq ../quast_results/latest/* 01_QUAST
rsync -avq ../F38-4-1.contigs.polished.fasta 02_Assembly
rsync -avq ../run_busco/*.tsv ../run_busco/*.txt 03_BUSCO
cd ..
