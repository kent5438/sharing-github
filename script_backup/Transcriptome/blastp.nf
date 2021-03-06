#!/usr/bin/env nextflow

params.nr = "/export/EC1680U/DataBase/Uniprot/Swiss-Prot/uniprot-reviewed_20190618.fasta"
params.nt = "/export/EC1680U/DataBase/blastdb/nt"
//params.gilist = "/export/EC1680U/DataBase/gilist/Chordata.prot.180802.gilist"
params.query = "4_Annotation/transdecoder/Trinity.95.fasta.transdecoder.pep"
params.out = "4_Annotation/blast_result/uniprot.blastp.outfmt6"
params.chunkSize = 2000

nrDB = file(params.nr)
//ntDB = file(params.nt)
//gilist = file(params.gilist)

Channel 
	.fromPath(params.query) 
	.splitFasta( by: params.chunkSize )
	.set { seq }



process blastp {

	cpus = 8
	executor = 'slurm'
	clusterOptions = '--job-name=blastp_nxf_kent'
//	clusterOptions = '--job-name=blastp_nxf_kent --partition=wgs'

    input:
    file 'seq.fa' from seq

    output:
    file 'out' into blast_result

    """
    /export/EC1680U/kentchen/miniconda3/envs/illumina-transcriptome/bin/blastp \
		-num_threads 8 -max_hsps 1 -max_target_seqs 1 -db $nrDB -query seq.fa \
		-outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle' \
		-evalue 1e-5 > out
    """
}

blast_result
	.collectFile(name: params.out)
	.subscribe { merged_file -> merged_file.copyTo(params.out) }

