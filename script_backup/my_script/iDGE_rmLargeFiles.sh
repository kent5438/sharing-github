#! /bin/bash

locate=`pwd | grep kentchen | grep "iDGE" | wc -l`
if [ ${locate} != 1 ]; then
	echo "### ERROR: You are not located at '/export/EC1680U/kentchen/iDGE'"
	exit 1;
fi

find . -type d -name work | xargs rm -rf
find . -type f -name '*.transcript.bam' | xargs rm -f
