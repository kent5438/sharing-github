#!/usr/bin/env nextflow

// located at 4_Annotation/transdecoder
params.db = "/export/EC1680U/DataBase/Trinotate/Pfam-A.hmm"
params.query = "4_Annotation/transdecoder/Trinity.95.fasta.transdecoder.pep"
params.out = "4_Annotation/features/PFAM.out"
params.chunkSize = 2000

DB = file(params.db)

Channel 
	.fromPath(params.query) 
	.splitFasta( by: params.chunkSize )
	.set { seq }



process hmmscan {

    clusterOptions = '--job-name=hmmscan_nxf_kent'

    input:
    file 'seq.fa' from seq

    output:
    file 'out' into hmmscan_result

    """
    /export/EC1680U/kentchen/miniconda3/envs/illumina-transcriptome/bin/hmmscan -E 0.01 \
		--cpu 8 --domtblout out $DB seq.fa \
    """
}

hmmscan_result
	.collectFile(name: params.out)
	.subscribe { merged_file -> merged_file.copyTo(params.out) }

