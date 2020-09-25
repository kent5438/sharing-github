#! /bin/bash

mkdir genome/STAR

# build index
STAR --runThreadN 32 --runMode genomeGenerate --genomeFastaFiles genome/mm10.fa --genomeDir genome/STAR

# mapping (use for loop)
STAR --runThreadN 32 --genomeDir genome/STAR --sjdbGTFfile genome/refGene.gtf --readFilesIn reads/C1_R1_001.clean.fastq.gz reads/C1_R2_001.clean.fastq.gz --quantMode TranscriptomeSAM --readFilesCommand zcat

mv Aligned.toTranscriptome.out.bam C1.bam

STAR --runThreadN 32 --genomeDir genome/STAR --sjdbGTFfile genome/refGene.gtf --readFilesIn reads/C2_R1_001.clean.fastq.gz reads/C2_R2_001.clean.fastq.gz --quantMode TranscriptomeSAM --readFilesCommand zcat

mv Aligned.toTranscriptome.out.bam C2.bam

STAR --runThreadN 32 --genomeDir genome/STAR --sjdbGTFfile genome/refGene.gtf --readFilesIn reads/S1_R1_001.clean.fastq.gz reads/S1_R2_001.clean.fastq.gz --quantMode TranscriptomeSAM --readFilesCommand zcat

mv Aligned.toTranscriptome.out.bam S1.bam

STAR --runThreadN 32 --genomeDir genome/STAR --sjdbGTFfile genome/refGene.gtf --readFilesIn reads/S2_R1_001.clean.fastq.gz reads/S2_R2_001.clean.fastq.gz --quantMode TranscriptomeSAM --readFilesCommand zcat

mv Aligned.toTranscriptome.out.bam S2.bam

### for loop example
#for i in `ls reads | cut -d'_' -f1 | sort -u`; do STAR --runThreadN 32 --genomeDir genome/STAR --sjdbGTFfile genome/CM.gtf --readFilesIn reads/${i}_R1_001.clean.fastq.gz reads/${i}_R2_001.clean.fastq.gz --quantMode TranscriptomeSAM --readFilesCommand zcat; rm -f Aligned.out.sam; mv Aligned.toTranscriptome.out.bam ${i}.bam; done
###

# merge BAM
picard MergeSamFiles I=C1.bam I=C2.bam I=S1.bam I=S2.bam O=all.sorted.bam SO=coordinate

# trinity genome-guided assembly
Trinity --genome_guided_bam bam/all.sorted.bam --genome_guided_max_intron 10000 --max_memory 200G --CPU 32


