#!/usr/bin/env nextflow

// located at 4_Annotation/transdecoder
params.db = "/export/EC1680U/DataBase/Trinotate/Pfam-A.hmm"
params.query = "Trinity.95.fasta.transdecoder_dir/longest_orfs.pep"
params.out = "pfam.domtblout"
params.chunkSize = 2000

DB = file(params.db)

Channel 
	.fromPath(params.query) 
	.splitFasta( by: params.chunkSize )
	.set { seq }



process hmmscan {

	cpus = 8
	executor = 'slurm'
//	clusterOptions = '--job-name=hmmscan_nxf_kent --partition=wgs'
    clusterOptions = '--job-name=hmmscan_nxf_kent'

    input:
    file 'seq.fa' from seq

    output:
    file 'out' into hmmscan_result

    """
    /export/EC1680U/kentchen/miniconda3/envs/illumina-transcriptome/bin/hmmscan \
		--cpu 8 --domtblout out $DB seq.fa \
    """
}

hmmscan_result
	.collectFile(name: params.out)
	.subscribe { merged_file -> merged_file.copyTo(params.out) }

