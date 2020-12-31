#!/usr/bin/env nextflow

// located at 4_Annotation/transdecoder
params.nr = "/export/EC1680U/DataBase/Uniprot/Swiss-Prot/uniprot-reviewed_20190618.fasta"
params.query = "Trinity.95.fasta.transdecoder_dir/longest_orfs.pep"
params.out = "blastp.outfmt6"
params.chunkSize = 2000

nrDB = file(params.nr)

Channel 
	.fromPath(params.query) 
	.splitFasta( by: params.chunkSize )
	.set { seq }



process blastp {

	clusterOptions = '--job-name=blastp_nxf_kent'

    input:
    file 'seq.fa' from seq

    output:
    file 'out' into blast_result

    """
    /export/EC1680U/kentchen/miniconda3/envs/illumina-transcriptome/bin/blastp \
		-num_threads 8 -max_target_seqs 1 -db $nrDB -query seq.fa \
		-outfmt 6 \
		-evalue 1e-5 > out
    """
}

blast_result
	.collectFile(name: params.out)
	.subscribe { merged_file -> merged_file.copyTo(params.out) }

