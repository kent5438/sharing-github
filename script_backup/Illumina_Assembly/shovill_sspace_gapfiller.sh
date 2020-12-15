#! /bin/bash
### Goals: (contigs.fasta) ---> (contigs.polished.fasta)

#now="[`date +"%F %T"`]"

prefix=$1
if [ -z ${prefix} ]; then
	echo "### ERROR: Give the prefix of sample name!!!
usage: $0 <prefix>";
	exit 1;
fi

# Main Pipeline

echo "Activate shovill environment..." | ts '[%Y-%m-%d %H:%M:%S]'
. activate shovill

echo "Check shovill contigs.fa existed or not..." | ts '[%Y-%m-%d %H:%M:%S]'
if [ ! -f "shovill/contigs.fa" ]; then
	shovill --tmpdir /home/kentchen/tmp --R1 reads/${prefix}_R1.clean.fastq.gz --R2 reads/${prefix}_R2.clean.fastq.gz --outdir shovill
fi

cd shovill

echo "lib01 bwa ../reads/${prefix}_R1.clean.fastq.gz ../reads/${prefix}_R2.clean.fastq.gz 414 0.5 FR" > lib.txt

echo "Scaffolding..." | ts '[%Y-%m-%d %H:%M:%S]'
if [ ! -f "standard_output/standard_output.final.scaffolds.fasta" ]; then
	perl /export/EC1680U/kentchen/opt/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl -l lib.txt -s contigs.fa -T 32
fi

echo "Gapfilling..." | ts '[%Y-%m-%d %H:%M:%S]'
if [ ! -f "standard_output/standard_output.gapfilled.final.fa" ]; then
	perl /export/EC1680U/kentchen/opt/GapFiller_v1-10_linux-x86_64/GapFiller.pl -l lib.txt -s standard_output/standard_output.final.scaffolds.fasta -T 32
fi

cd ..

echo "Finish ~ <(_ _)>" | ts '[%Y-%m-%d %H:%M:%S]'
