#! /bin/bash

longest_orf='Trinity.95.fasta.transdecoder_dir/longest_orfs.pep'
pred_orf='Trinity.95.fasta.transdecoder.pep'

if [ ! -e "4_Annotation/transdecoder" ]; then
	mkdir -p 4_Annotation/transdecoder
fi

if [ ! -e "4_Annotation/blast_result" ]; then
	mkdir -p 4_Annotation/blast_result
fi

if [ ! -e "4_Annotation/features" ]; then
	mkdir -p 4_Annotation/features
fi

cd 4_Annotation/transdecoder/

if [ ! -e "$longest_orf" ]; then
	TransDecoder.LongOrfs -t ../../1_AssemblyStats/Trinity.95.fasta
fi

nextflow run /export/EC1680U/perl/bin/Trinity/blastp_transdecoder.nf -c /export/EC1680U/perl/bin/Trinity/nextflow.config -resume
nextflow run /export/EC1680U/perl/bin/Trinity/hmmscan_transdecoder.nf -c /export/EC1680U/perl/bin/Trinity/nextflow.config -resume


if [ ! -e "$pred_orf" ]; then
	TransDecoder.Predict -t ../../1_AssemblyStats/Trinity.95.fasta \
		--retain_pfam_hits pfam.domtblout \
		--retain_blastp_hits blastp.outfmt6
fi

rm -f pfam.domtblout
rm -f blastp.outfmt6

cd ../..
