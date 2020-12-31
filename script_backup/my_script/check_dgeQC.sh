#! /bin/bash
printf "[Sequencing Total Bases]\n"
grep 'total bases' results/4.stats_out/*.fastq-stats.txt | cut -d'/' -f3

printf "\n"

printf "[# of Paired reads]\n"
grep -P 'reads\t' results/4.stats_out/*.fastq-stats.txt | cut -d'/' -f3

printf "\n"

printf "[Alignment Summary]\n"
grep -i 'overall' results/1.CalculateExpression/.*.command.err | cut -d'/' -f3 | sed 's/^\.//'
