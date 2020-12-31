#!/usr/bin/env nextflow

params.nr = "/export/EC1680U/DataBase/Uniprot/Swiss-Prot/uniprot-reviewed_20190618.fasta"
params.nt = "/export/EC1680U/DataBase/blastdb/nt"
//params.gilist = "/export/EC1680U/DataBase/gilist/Chordata.prot.180802.gilist"
params.query = "1_AssemblyStats/Trinity.95.fasta"
params.out = "4_Annotation/blast_result/uniprot.blastx.outfmt6"
params.chunkSize = 2000

nrDB = file(params.nr)
//ntDB = file(params.nt)
//gilist = file(params.gilist)

Channel 
	.fromPath(params.query) 
	.splitFasta( by: params.chunkSize )
	.set { seq }


process blastx {

	clusterOptions = '--job-name=blastx_nxf_kent'

    input:
    file 'seq.fa' from seq
//    file 'gilist' from gilist

    output:
    file 'out' into blast_result

    """
    /export/EC1680U/kentchen/miniconda3/envs/illumina-transcriptome/bin/blastx \
		-num_threads 8 -max_hsps 1 -max_target_seqs 1 -db $nrDB -query seq.fa \
		-outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle' > out
#		-outfmt '6 qseqid sseqid pident evalue stitle' > out
    """
}

blast_result
	.collectFile(name: params.out)
	.subscribe { merged_file -> merged_file.copyTo(params.out) }

