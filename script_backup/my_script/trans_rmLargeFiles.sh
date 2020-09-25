#! /bin/bash

locate=`pwd | grep kentchen | grep "Transcriptome" | wc -l`
if [ ${locate} != 1 ]; then
    echo "### ERROR: You are NOT located at '/export/EC1680U/kentchen/Transcriptome'"
fi

find . -name '*.minC1.gff' | xargs rm -rf
find . -name 'all.sorted.bam.*' | xargs rm -rf
find . -name '*.sam' | xargs rm -rf
find . -name 'chrysalis' | xargs rm -rf
find . -name 'insilico_read_normalization' | xargs rm -rf
find . -name 'jellyfish.kmers.25.asm.fa' | xargs rm -rf
find . -name 'read_partitions' | xargs rm -rf
find . -name 'both.fa' | xargs rm -rf
find . -name 'bowtie2.bam' | xargs rm -rf
find . -name 'bowtie2.bam.for_rsem.bam' | xargs rm -rf
find . -name 'all.sorted.bam' | xargs rm -rf
