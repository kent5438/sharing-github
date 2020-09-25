#! /bin/sh

for i in *.fastq.gz; do mv $i `echo $i | awk -F'_' '{print $1 "_" $4}'`; done
