#! /bin/bash

locate=`pwd | grep kentchen | grep "16S" | wc -l`
if [ ${locate} != 1 ]; then
	echo "### ERROR: You are NOT located at '/export/EC1680U/kentchen/iDGE'"
fi

find . -type d -name work | xargs sudo rm -rf
find . -type f -name '*.unique.align' | xargs rm -f
find . -type f -name '*.unique.dist' | xargs rm -f
