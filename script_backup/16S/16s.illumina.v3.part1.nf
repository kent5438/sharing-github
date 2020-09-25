#!/usr/bin/env nextflow

/********************************************************************************
 * author: 劉丞原 (daniel.liu)
 * Release date: 2019/01/09
 * Modified date:
 * nextflow run 16s.illumina.v2.nf --rm -c 16.config -with-report report.html -with-dag 16s.svg
 * 16S 的 nextflow 加上了不同運算環境的設定檔
 * 執行方式如下：
 * 1. configure file為 /export/EC1680U/perl/bin/16S/nextflow.manually.config
 * 2. profile： 目前可用的為gbc, singularity, docker, standard
 * gbc: 基米slurm，queue為all
 * standard: local執行
 * 舉例 -profile standard,singularity    在local執行並使用singularity image
 *     -profile gbc,singularity             透過slurm分配工作，使用singularity image
 * nextflow run /export/EC1680U/perl/bin/16S/16s.illumina.v3.nf -c /export/EC1680U/perl/bin/16S/nextflow.manually.config -profile standard,docker -with-report report.html
********************************************************************************/
version = '3.0'

params.file = "Mapping.txt"
params.group = 'Description'
params.lefse_group = 'Description'
if (params.file) {
  mapping_file = file(params.file)
} else {
  mapping_file = file("Mapping.txt")
}

if (params.outdir){
  directoryMap = defineDirectoryMap(params.outdir)
} else{
  exit 1, 'You should specifiy the output directory by using --outdir'
}

sample_reads = Channel
  .from(mapping_file.readLines())
  .map { line ->
    m = line =~ /(.*)\t(.*)\t(.*)\t(.*)\t(.*)/
    id = m[0][1]
    id = id.replaceAll(/-/, "_")
    if (id=~/#SampleID/) {
    } else {
      desc = m[0][4]
      filepath = m[0][5]
      m2 = filepath =~/(.*),(.*)/
      R1 = m2[0][1]
      R2 = m2[0][2]
      [desc, id, R1, R2]
    }
   }
(sample_reads_qc, sample_reads_trim) = sample_reads.into(2)
process QC_Raw_Data {
  publishDir directoryMap.qc_rawData, mode: 'copy', overwrite: true
  tag "$sampleid"

  input:
  set val(group), val(sampleid), val(r1), val(r2) from sample_reads_qc

  output:
  set val(group), val(sampleid), file('*_fastqc.{zip,html}') into fastqc_rawdata_output

  script:
  """
  fastqc -o . ${r1}
  fastqc -o . ${r2}
  """
}
process MQc_Raw_Data {
  publishDir directoryMap.mqc_rawdata, mode: 'copy', overwrite: true
  tag "MQcRawData - All"

  input:
  file ('B.1_QC/fastqc-rawdata/*') from fastqc_rawdata_output.flatten().toList()

  output:
  set file('multiqc_report.html'), file('multiqc_data') into mqc_rawdata_output

  script:
  """
  multiqc -f -q -c /export/EC1680U/perl/bin/16S/MiSeq-16S_multiqc_config.yaml -o . .
  """
}
process Trimming_Raw_Data {
  publishDir directoryMap.trimming_rawdata, mode: 'copy', overwrite: true
  tag "$sampleid"

  input:
  set val(group), val(sampleid), val(r1), val(r2) from sample_reads_trim

  output:
  set val(group), val(sampleid), file(trim_r1), file(trim_r2) into trim_rawdata_output

  script:
  trim_r1 = sampleid + ".R1.clean.fastq"
  trim_r2 = sampleid + ".R2.clean.fastq"

  """
  # CCTACGGGNGGCWGCAG Bakt_314F
  # GACTACHVGGGTATCTAATCC Bakt_805R
  cutadapt -m ${params.read_length} --discard-untrimmed -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -o ${trim_r1} -p ${trim_r2} ${r1} ${r2}
  """
}
(trim_rawdata_output_qc, trim_rawdata_output_join) = trim_rawdata_output.into(2)
process Qc_Trimed_Data {
  publishDir directoryMap.qc_trimeddata, mode: 'copy', overwrite: true
  tag "$sampleid"

  input:
  set val(group), val(sampleid), file(r1), file(r2) from trim_rawdata_output_qc

  output:
  set val(group), val(sampleid), file('*_fastqc.{zip,html}') into fastqc_trimdedata_output

  script:
  """
  fastqc -o . ${r1}
  fastqc -o . ${r2}
  """
}
process MQc_Trimed_Data {
  publishDir directoryMap.mqc_trimeddata, mode: 'copy', overwrite: true
  tag "MQcTrimedData - All"

  input:
  file ('B.1_QC/fastqc-trimmed/*') from fastqc_trimdedata_output.flatten().toList()

  output:
  set file('multiqc_report.html'), file('multiqc_data') into mqc_trimdedata_output

  script:
  """
  multiqc -f -q -c /export/EC1680U/perl/bin/16S/MiSeq-16S_multiqc_config.yaml -o . .
  """
}
process Join_Fastq {
  publishDir directoryMap.join_fastq, mode: 'copy', overwrite: true
  tag "$sampleid"
  cpus 20

  input:
  set val(group), val(sampleid), file(trim_r1), file(trim_r2) from trim_rawdata_output_join

  output:
  set val(group), val(sampleid), file("${sampleid}.extendedFrags.fastq") into join_fastq_output

  script:
  """
  flash ${trim_r1} ${trim_r2} -p ${params.phred_offset} -m ${params.min_overlap} -M ${params.max_overlap} -x ${params.mismatch_ratio} -d . -o ${sampleid} -t ${task.cpus} > ${sampleid}.log 2>&1;
  """
}
process Dem_Qua_Filter {
  publishDir directoryMap.dem_qua_fil, mode: 'copy', overwrite: true
  tag "$sampleid"

  input:
  set val(group), val(sampleid), file(fastq) from join_fastq_output

  output:
  set val(group), val(sampleid), file("${sampleid}.fna") into demultiplex_quality_filter_output

  script:
  output = sampleid + '.fna'
  """
  perl /export/EC1680U/perl/bin/16S/mergeSeq.nextflow.pl ${mapping_file} .
  mv seqs.fna ${output}
  """
}
(demultiplex_quality_filter_output1, demultiplex_quality_filter_output2) = demultiplex_quality_filter_output.into(2)
process Rem_Chimeras {
  publishDir directoryMap.rem_chimeras, mode: 'copy', overwrite: true
  tag "$sampleid"
  cpus 20

  input:
  set val(group), val(sampleid), file(fna) from demultiplex_quality_filter_output1

  output:
  set val(group), val(sampleid), file("${sampleid}.trimmed.fa") into removing_chimeras_output

  script:

  """
  export TERM=xterm
  /app_tools/mothur/mothur "#chimera.uchime(processors=${task.cpus}, fasta=${fna}, reference=${params.chimeras_ref})" > /dev/null 2>&1
  removing_chimera.pl '${fna}' '${sampleid}.ref.uchime.chimeras' '${sampleid}.trimmed.fa'
  """
}
process Otu_Preprocess {
  publishDir directoryMap.otu_preprocess, mode: 'copy', overwrite: true
  tag "$sampleid"
  cpus 20

  input:
  set val(group), val(sampleid), file(fa) from removing_chimeras_output

  output:
  set val(group), val(sampleid), file("${sampleid}.trimmed.unique.fa"), file("${sampleid}.trimmed.count_table"), file("${sampleid}.trimmed.unique.align"),  file("${sampleid}.trimmed.unique.silva.wang.taxonomy") into otu_output

  script:
  """
  export TERM=xterm
  mothur_make_group.pl ${fa}
  /app_tools/mothur/mothur "#unique.seqs(fasta=${fa})" > /dev/null # 輸出 XX.trimmed.unique.fa & XX.trimmed.names
  /app_tools/mothur/mothur "#count.seqs(name=${sampleid}.trimmed.names, group=group.txt)" > /dev/null # 輸出 XX.trimmed.count_table
  /app_tools/mothur/mothur "#classify.seqs(processors=${task.cpus}, fasta=${sampleid}.trimmed.unique.fa, name=${sampleid}.trimmed.names, template=${params.silva_otu_ref_fasta}, taxonomy=${params.silva_otu_ref_taxonomy}, cutoff=0, method=wang)"  # 輸出 XX.trimmed.unique.raw_taxonomy.wang.taxonomy  & XX.trimmed.unique.raw_taxonomy.wang.tax.summary
  /app_tools/mothur/mothur "#remove.lineage(fasta=${sampleid}.trimmed.unique.fa, count=${sampleid}.trimmed.count_table, taxonomy=${sampleid}.trimmed.unique.raw_taxonomy.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-mitochondria-chloroplast)" > /dev/null # 輸出 XX.trimmed.unique.raw_taxonomy.wang.pick.taxonomy & XX.trimmed.unique.pick.fa & XX.trimmed.pick.count_table
  mv ${sampleid}.trimmed.unique.raw_taxonomy.wang.pick.taxonomy ${sampleid}.trimmed.unique.silva.wang.taxonomy
  mv ${sampleid}.trimmed.unique.pick.fa ${sampleid}.trimmed.unique.fa
  mv ${sampleid}.trimmed.pick.count_table ${sampleid}.trimmed.count_table
  /app_tools/mothur/mothur "#align.seqs(processors=${task.cpus}, candidate=${sampleid}.trimmed.unique.fa, template=${params.silva_otu_ref_aligned_fasta})" > /dev/null # 輸出 XX.trimmed.unique.align & XX.trimmed.unique.align.report && XX.trimmed.unique.flip.accnos
  """
}
process Otu_Similarity97 {
  publishDir directoryMap.otu_preprocess, mode: 'copy', overwrite: true
  tag "$sampleid"
  cpus 20

  input:
  set val(group), val(sampleid), file(fa), file(count_table), file(align), file(taxonomy) from otu_output

  output:
  set val(group), val(sampleid), file("97/${sampleid}.trimmed.unique.opti_mcc.list"), file("97/${sampleid}.trimmed.unique.opti_mcc.0.03.cons.taxonomy"),  file("97/${sampleid}.otu_table.taxID.filter.dense.dedup.biom"), file("97/${sampleid}.otu_table.taxID.filter.dense.dedup.txt"), file('97/*') into otu_output_similarity97

  script:
  """
  export TERM=xterm
  mkdir 97
  cp ${fa} 97
  cp ${count_table} 97
  cp ${align} 97
  cp ${taxonomy} 97
  cd 97/
  /app_tools/mothur/mothur "#dist.seqs(processors=${task.cpus}, fasta=${align}, countends=F, cutoff=0.03)" > /dev/null
  /app_tools/mothur/mothur "#cluster(column=${sampleid}.trimmed.unique.dist, count=${count_table}, method=opti, cutoff=0.03)" > /dev/null
  /app_tools/mothur/mothur "#make.shared(list=${sampleid}.trimmed.unique.opti_mcc.list, count=${count_table})" > /dev/null
  /app_tools/mothur/mothur "#classify.otu(list=${sampleid}.trimmed.unique.opti_mcc.list, count=${count_table}, cutoff=80, label=0.03, taxonomy=$taxonomy)'" > /dev/null
  #grep 'unknown;' ${sampleid}.trimmed.unique.opti_mcc.0.03.cons.taxonomy | cut -f1 > unknown_ids.txt
  ## OTUs ref ID, used for comparing with OTU sequence
  /app_tools/mothur/mothur "#make.biom(shared=${sampleid}.trimmed.unique.opti_mcc.shared, constaxonomy=${sampleid}.trimmed.unique.opti_mcc.0.03.cons.taxonomy, matrixtype=dense)" > /dev/null
  mv ${sampleid}.trimmed.unique.opti_mcc.0.03.biom ${sampleid}.otu_table.otuID.biom
  #filter_otus_from_otu_table.py -i ${sampleid}.otu_table.otuID.biom -e unknown_ids.txt -o ${sampleid}.otu_table.otuID.filter.biom --min_count_fraction 0.00005 # ref: https://goo.gl/EMwR9N
  filter_otus_from_otu_table.py -i ${sampleid}.otu_table.otuID.biom -o ${sampleid}.otu_table.otuID.filter.biom --min_count_fraction 0.00005 # ref: https://goo.gl/EMwR9N
  biom convert -i ${sampleid}.otu_table.otuID.filter.biom -o ${sampleid}.otu_table.otuID.filter.txt --to-tsv --header-key taxonomy --table-type="OTU table"

  ## Generate OTUs sequence and taxonomy output
  /app_tools/mothur/mothur "#get.oturep(column= ${sampleid}.trimmed.unique.dist, list= ${sampleid}.trimmed.unique.opti_mcc.list, count=${count_table}, fasta=${fa})" > /dev/null

  ## Convert otuID.filter to taxID.filter
  vsearch --usearch_global  ${sampleid}.trimmed.unique.opti_mcc.0.03.rep.fasta --db ${params.silva_otu_ref_fasta} --id 0.1 --iddef 1 --blast6out results.blast6 --maxaccepts 1 --threads 32


  ### Contributed by JSYF ###
  perl /export/EC1680U/perl/bin/16S/biomProc.pl sparseToDense \
  	${sampleid}.otu_table.otuID.filter.biom \
  	${sampleid}.otu_table.otuID.filter.dense.biom

  perl /export/EC1680U/perl/bin/16S/biomProc.pl idChangeUsingFasta \
  	${sampleid}.otu_table.otuID.filter.dense.biom \
  	${sampleid}.otu_table.taxID.filter.dense.biom \
  	${sampleid}.trimmed.unique.opti_mcc.0.03.rep.fasta,results.blast6

  perl /export/EC1680U/perl/bin/16S/biomProc.pl uniqById \
  	${sampleid}.otu_table.taxID.filter.dense.biom \
  	${sampleid}.otu_table.taxID.filter.dense.dedup.biom
  ########################

  biom convert -i ${sampleid}.otu_table.taxID.filter.dense.dedup.biom -o ${sampleid}.otu_table.taxID.filter.dense.dedup.txt --to-tsv --header-key taxonomy --table-type="OTU table"
  rm ${fa}
  rm ${count_table}
  rm ${align}
  rm ${taxonomy}
  rm results.blast6
  rm mothur.*.logfile
  cd ../
  """
}
def defineDirectoryMap(outdir) {
  return [
    'root_dir'         : "$outdir/",
    'qc_rawData'       : "$outdir/B.1_QC/fastqc-rawdata",
    'mqc_rawdata'      : "$outdir/B.1_QC/multiqc-rawdata",
    'trimming_rawdata' : "$outdir/B.1_QC/trimming",
    'qc_trimeddata'    : "$outdir/B.1_QC/fastqc-trimmed",
    'mqc_trimeddata'   : "$outdir/B.1_QC/multiqc-trimmed",
    'join_fastq'       : "$outdir/B.2_join_fastq",
    'dem_qua_fil'      : "$outdir/B.3_demultiplex_quilty_filter",
    'rem_chimeras'     : "$outdir/B.4_removing_chimeras",
    'otu_preprocess'   : "$outdir/B.5_OTU_cluster",
    'stat_analysis'    : "$outdir/B_statistic_analysis",
    'alpha_diversity'  : "$outdir/C_alpha_diversity",
    'rank_abundance'   : "$outdir/C.4_rank_abundance",
    'alpha_indices_table' : "$outdir/C_alpha_indices_table",
    'taxa_summary'     : "$outdir/D.1_taxa_summary",
    'new_taxa_summary' : "$outdir/D.1_new_taxa_summary",
    'translated'       : "$outdir/D.1_translated",
    'heatmap'          : "$outdir/D.2_heatmap",
    'phylogenetic_tree': "$outdir/D.3_phylogenetic_tree",
    'metastats'        : "$outdir/D.4_metastats",
    'beta_diversity'   : "$outdir/E.1_PCA_PCoA",
    'upgma'            : "$outdir/E.2_UPGMA",
    'lefse'            : "$outdir/E.3_LEfSe",
    'nmds'             : "$outdir/F.1_NMDS",
    'anosim'           : "$outdir/F.3_Anosim",
    'mrpp'             : "$outdir/F.4_MRPP",
    'dca_rda_cca'      : "$outdir/F_3A",
    'tax4fun'          : "$outdir/F.5_Tax4Fun",
  ]
}

workflow.onComplete {
    log.info "Work Dir    : $workflow.workDir"
    log.info "Completed at: $workflow.complete"
    log.info "Duration    : $workflow.duration"
    log.info "Success     : $workflow.success"
    log.info "Exit status : $workflow.exitStatus"
    log.info "Error report: " + (workflow.errorReport ?: '-')
    log.info "Command line: $workflow.commandLine"
}

def returnFile(it) {
  // return file if it exists
  final f = file(it)
  if (!f.exists()) {
    exit 1, "Missing file "
  }
  return f
}
workflow.onError {
  log.info "Workflow execution stopped with the following message: " + workflow.errorMessage
}
