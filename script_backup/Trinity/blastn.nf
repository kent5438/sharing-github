#!/usr/bin/env nextflow

params.nt = "/export/EC1680U/DataBase/blastdb/nt"
params.gilist = "/export/EC1680U/DataBase/gilist/Tetranychidae.nucl.181005.gilist"
params.query = "1_AssemblyStats/Trinity.95.fasta"
params.out = "uniprot.blastn.outfmt6"
params.chunkSize = 2000

//nrDB = file(params.nr)
ntDB = file(params.nt)
gilist = file(params.gilist)

Channel 
	.fromPath(params.query) 
	.splitFasta( by: params.chunkSize )
	.set { seq }



process blastn {

	cpus = 8
	executor = 'slurm'
	clusterOptions = '--job-name=blastn_nxf_kent'

    input:
    file 'seq.fa' from seq
    file 'gilist' from gilist

    output:
    file 'out' into blast_result

    """
    /export/EC1680U/kentchen/miniconda3/envs/illumina-transcriptome/bin/blastn \
        -num_threads 8 -max_target_seqs 1 -db $ntDB -gilist $gilist -query seq.fa \
        -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle' \
        -evalue 1e-6 > out
    """

}

blast_result
	.collectFile(name: params.out)
	.subscribe { merged_file -> merged_file.copyTo(params.out) }

