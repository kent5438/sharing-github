#!/usr/bin/env nextflow

//params.long = "PS18014-PS19017_merged.fasta"
//params.short = "HS19036.interleaved.fastq"
//params.out = "cumcumis_cor"

long_file = file(params.long_file)
short_file = file(params.short_file)
out = file(params.out)

Channel 
	.fromPath( long_file )
	.splitFasta( by: 2000 )
	.set { chunk }


process lordec {

	cpus = 32
	executor = 'slurm'
	clusterOptions = '--job-name=lordec_nxf_kent'

    input:
    file 'chunks.fa' from chunk
    file 'short.fq' from short_file

    output:
    file 'outs' into result

    """
	lordec-correct -i chunks.fa -2 short.fq -o outs -k 19 -s 3 -T 32
	rm -f short.fq_k19_s3.h5
    """

}

result
	.collectFile(name: out)
	.subscribe { merged_file -> merged_file.copyTo(out) }

