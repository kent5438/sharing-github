#! /bin/bash

samtools faidx $1
cut -f2 ${1}.fai | Rscript -e 'data <- as.numeric (readLines ("stdin")); summary(data); hist(data)'
