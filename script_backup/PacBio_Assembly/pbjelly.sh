#! /bin/bash

. activate pbjelly

source /export/EC1680U/kentchen/opt/PBSuite_15.8.24/setup.sh

ID=$1

if [ ! -d "pbjelly_out" ]; then
	mkdir pbjelly_out
fi

if [ ! -f "subreads.fasta" ]; then
	if [ -f "../../subreads/${ID}.subreads.fastq" ]; then
		seqtk seq -a ../../subreads/${ID}.subreads.fastq > subreads.fasta
	else
		ln -s ../../subreads/${ID}.subreads.fasta subreads.fasta
	fi
fi

cp /export/EC1680U/kentchen/opt/PBSuite_15.8.24/Protocol.xml .
sed -i "s#referenceFasta#${PWD}/scaffolds.fasta#" Protocol.xml
sed -i "s#outputPath#${PWD}/pbjelly_out#" Protocol.xml
sed -i "s#readdirectory#${PWD}#" Protocol.xml

if [ ! -f "subreads.qual" ]; then
	fakeQuals.py subreads.fasta subreads.qual
fi

Jelly.py setup Protocol.xml
Jelly.py mapping Protocol.xml
Jelly.py support Protocol.xml
Jelly.py extraction Protocol.xml
Jelly.py assembly Protocol.xml
Jelly.py output Protocol.xml

