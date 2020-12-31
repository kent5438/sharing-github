#! /bin/bash

locate=`pwd | grep kentchen | grep "Pacbio" | wc -l`
if [ ${locate} != 1 ]; then
	echo "### ERROR: You are not located at '/export/EC1680U/kentchen/Pacbio'"
	exit 1;
fi

find . -type f -name 'subreads.qual' | xargs rm -f
find . -type f -name '*.sorted.bam' | xargs rm -f
find . -type f -name '*.subreads.fasta' | xargs rm -f
