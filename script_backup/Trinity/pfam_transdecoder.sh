#! /bin/bash
hmmscan --cpu 32 --domtblout pfam.domtblout \
	/export/EC1680U/DataBase/Trinotate/Pfam-A.hmm \
	Trinity.95.fasta.transdecoder_dir/longest_orfs.pep
