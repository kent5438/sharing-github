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

params.file = "data/Mapping.txt"
params.group = 'Description'
params.lefse_group = 'Description'
if (params.file) {
  mapping_file = file(params.file)
} else {
  mapping_file = file("data/Mapping.txt")
}

if (params.outdir){
  directoryMap = defineDirectoryMap(params.outdir)
} else{
  exit 1, 'You should specifiy the output directory by using --outdir'
}

sample_reads = Channel
  .from(mapping_file.readLines())
  .map { line ->
    m = line =~ /^(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S.*)/
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
  # FWD: CCTACGGGNGGCWGCAG Bakt_314F
  # REV: GACTACHVGGGTATCTAATCC Bakt_805R
  # FWD.RC: CTGCWGCCNCCCGTAGG
  # REV.RC: GGATTAGATACCCBDGTAGTC
  cutadapt -m ${params.read_length} --discard-untrimmed -g CCTACGGGNGGCWGCAG -g GACTACHVGGGTATCTAATCC -G CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -o ${trim_r1} -p ${trim_r2} ${r1} ${r2}
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
#  /app_tools/mothur/mothur "#chimera.uchime(processors=${task.cpus}, fasta=${fna}, reference=${params.chimeras_ref})" > /dev/null 2>&1
  /export/EC1680U/kentchen/opt/usearch11.0.667_i86linux32 -uchime2_ref ${fna} -db ${params.chimeras_ref} -strand plus -threads ${task.cpus} -mode high_confidence -chimeras ${sampleid}.ref.uchime.chimeras -notmatched ${sampleid}.trimmed.fa
#  removing_chimera.pl '${fna}' '${sampleid}.ref.uchime.chimeras' '${sampleid}.trimmed.fa'
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
