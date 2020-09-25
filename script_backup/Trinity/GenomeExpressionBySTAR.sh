#! /bin/bash
#################################################################################
# Designed for Trinity analysis													#
# This pipeline is used for mapping reads to genome.							#
# It is not a necessary work if customer is satisified with transcript level.	#
# The downstream variant calling is concated with GenomeVariantBySTAR.sh		#
#################################################################################

DIR=`pwd`
REF=$1
SAMPLE=$2
GTF=${REF}.gtf
FNA=${REF}.fna

if [ "$1" == "-h" ]; then
	echo "Usage: Please follow the steps:
1. Please locate at project root directory...
2. Make sure genome/, genome/RSEM_ref/, genome/STAR/ folders exist
3. Make sure *.fna & *.gtf exist in genome/
4. Make sure 2_Quantitation/ exist
5. run % GenomeExpressionBySTAR.sh $REF $SAMPLE"
	exit 0
fi

# Convert gff to gtf
if [ ! -e "genome/${REF}.gtf" ] && [ -e "genome/${REF}.gff" ]; then
	/export/EC1680U/software/cufflinks-2.2.1.Linux_x86_64/gffread ${REF}.gff \
		-T -o ${REF}.gtf
fi

# Check ${REF}.cor.gtf
if [ ! -e "genome/${REF}.cor.gtf" ]; then
	perl -ne 'chomp; @a=split/\t/; %h=split(/ /,$a[8]); $a[8]=join(" ",("gene_id",$h{"gene_name"},"transcript_id",$h{"transcript_id"})); print join("\t",@a),"\n";' genome/${REF}.gtf > genome/${REF}.cor.gtf
fi


mkdir -p genome/RSEM_ref
mkdir -p genome/STAR

# build RSEM index
if [ ! -e "${DIR}/genome/RSEM_ref/RSEM_ref.n2g.idx.fa" ]; then
/export/EC1680U/software/anaconda2/bin/rsem-prepare-reference --gtf genome/${REF}.cor.gtf genome/${REF}.fna genome/RSEM_ref/RSEM_ref
fi

# build STAR index
if [ ! -e "${DIR}/genome/STAR/SAindex" ]; then
STAR --runThreadN 32 --runMode genomeGenerate --genomeFastaFiles genome/${REF}.fna --sjdbGTFfile genome/${REF}.cor.gtf --genomeDir genome/STAR/
fi

mkdir -p 2_Quantitation/${SAMPLE}/STAR
cd 2_Quantitation/${SAMPLE}/STAR

# mapping reads to transcripts
if [ ! -e "${DIR}/2_Quantitation/${SAMPLE}/STAR/Aligned.toTranscriptome.out.bam" ]; then
STAR --runThreadN 32 --genomeDir ${DIR}/genome/STAR/ --sjdbGTFfile ${DIR}/genome/${REF}.cor.gtf --readFilesIn ${DIR}/reads/${SAMPLE}*R1*fastq ${DIR}/reads/${SAMPLE}*R2*fastq --quantMode TranscriptomeSAM 
#STAR --runThreadN 32 --genomeDir ${DIR}/genome/STAR/ --sjdbGTFfile ${DIR}/genome/${REF}.cor.gtf --readFilesIn ${DIR}/reads/${SAMPLE}*R1*fastq.gz ${DIR}/reads/${SAMPLE}*R2*fastq.gz --quantMode TranscriptomeSAM --readFilesCommand zcat
fi

# mapping transcript-bam to genome 
if [ ! -e "${DIR}/2_Quantitation/${SAMPLE}/STAR/${SAMPLE}.isoforms.results" ]; then
/export/EC1680U/software/anaconda2/bin/rsem-calculate-expression --bam -p 32 --paired-end Aligned.toTranscriptome.out.bam ${DIR}/genome/RSEM_ref/RSEM_ref ${SAMPLE} --output-genome-bam
fi

rm -rf _STARgenome *.sam

cd $PWD
