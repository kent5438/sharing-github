#!/usr/bin/env nextflow

/********************************************************************************
 * author: 劉丞原 (daniel.liu)
 * Release date: 2018/03/22
 * Modified date:
 * nextflow run test.nf --rm  -c 16.config
 * nextflow run test.nf --rm -file NS18014/Mutect2.json -c gatk-mutect2.config --outdir NS18014 -w work/gatk4 -with-dag gatk-mutect2.svg
 ********************************************************************************/
version = '1.0'

params.file = "data/Mapping.txt"
params.group = 'Description'
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

/*
# restart Otu_similarity97
list = Channel
  .from(mapping_file.readLines())
  .map { line ->
    m = line =~ /(.*)\t(.*)\t(.*)\t(.*)\t(.*)/
    id = m[0][1]
    id = id.replaceAll(/-/, "_")
    if (id=~/#SampleID/) {
    } else {
      desc = m[0][4]
      table = "/export/EC1680U/Daniel_16S/16s/OUTPUT/B.5_OTU_cluster/" + id + ".trimmed.count_table"
      align = "/export/EC1680U/Daniel_16S/16s/OUTPUT/B.5_OTU_cluster/" + id + ".trimmed.unique.align"
      taxon = "/export/EC1680U/Daniel_16S/16s/OUTPUT/B.5_OTU_cluster/" + id + ".trimmed.unique.silva.wang.taxonomy"
      [desc, id, table, align, taxon]
    }
  }

if (params.step == 'similarity97')  {
   otu_output = list.map {
      it-> [it[0], it[1], file(it[2]), file(it[3]), file(it[4])]
   }
}
# restart Otu_preprocess
list = Channel
  .from(mapping_file.readLines())
  .map { line ->
    m = line =~ /(.*)\t(.*)\t(.*)\t(.*)\t(.*)/
    id = m[0][1]
    id = id.replaceAll(/-/, "_")
    if (id=~/#SampleID/) {
    } else {
      desc = m[0][4]
      fa = "/export/EC1680U/Daniel_16S/16s/OUTPUT/B.4_removing_chimeras/" + id + ".trimmed.fa"
      [desc, id, fa]
    }
  }
if (params.step == 'Otu_preprocess')  {
   removing_chimeras_output = list.map {
      it-> [it[0], it[1], file(it[2])]
   }
}
*/
process Qc_Raw_Data {
  publishDir directoryMap.Qc_RawData, mode: 'copy', overwrite: true
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
  publishDir directoryMap.MQc_RawData, mode: 'copy', overwrite: true
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
  publishDir directoryMap.Trimming_RawData, mode: 'copy', overwrite: true
  tag "$sampleid"

  input:
  set val(group), val(sampleid), val(r1), val(r2) from sample_reads_trim

  output:
  set val(group), val(sampleid), file(trim_r1), file(trim_r2) into trim_rawdata_output

  script:
  trim_r1 = sampleid + ".R1.clean.fastq"
  trim_r2 = sampleid + ".R2.clean.fastq"

  """
  cutadapt -m ${params.read_length} --discard-untrimmed -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -o ${trim_r1} -p ${trim_r2} ${r1} ${r2} > /dev/null
  """
}
(trim_rawdata_output_qc, trim_rawdata_output_join) = trim_rawdata_output.into(2)
process Qc_Trimed_Data {
  publishDir directoryMap.Qc_TrimedData, mode: 'copy', overwrite: true
  tag "$sampleid"

  input:
  set val(group), val(sampleid), val(r1), val(r2) from trim_rawdata_output_qc

  output:
  set val(group), val(sampleid), file('*_fastqc.{zip,html}') into fastqc_trimdedata_output

  script:
  """
  fastqc -o . ${r1}
  fastqc -o . ${r2}
  """
}
process MQc_Trimed_Data {
  publishDir directoryMap.MQc_TrimedData, mode: 'copy', overwrite: true
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
  publishDir directoryMap.Join_Fastq, mode: 'copy', overwrite: true
  tag "$sampleid"
  cpus 32

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
  publishDir directoryMap.Dem_Qua_Fil, mode: 'copy', overwrite: true
  tag "$sampleid"

  input:
  set val(group), val(sampleid), file(fastq) from join_fastq_output

  output:
  set val(group), val(sampleid), file("${sampleid}.histograms.txt"), file("${sampleid}.fna") into demultiplex_quality_filter_output

  script:

  output1 = sampleid + '.histograms.txt'
  output2 = sampleid + '.fna'
  """
  split_libraries_fastq.py --phred_quality_threshold ${params.quality} --phred_offset ${params.phred_offset} --barcode_type not-barcoded -i ${fastq} --sample_ids ${sampleid} -o .
  mv histograms.txt ${output1}
  mv seqs.fna ${output2}
  """
}
(demultiplex_quality_filter_output1, demultiplex_quality_filter_output2) = demultiplex_quality_filter_output.into(2)
process Rem_Chimeras {
  publishDir directoryMap.Rem_Chimeras, mode: 'copy', overwrite: true
  tag "$sampleid"
  cpus 32

  input:
  set val(group), val(sampleid), file(txt), file(fna) from demultiplex_quality_filter_output1

  output:
  set val(group), val(sampleid), file("${sampleid}.trimmed.fa") into removing_chimeras_output

  script:
  """
  #export PATH=/app_tools/mothur/:$PATH
  export TERM=xterm
  /app_tools/mothur/mothur "#chimera.uchime(processors=${task.cpus}, fasta=${fna}, reference=${params.chimeras_ref})" > /dev/null 2>&1
  removing_chimera.pl '${fna}' '${sampleid}.ref.uchime.chimeras' '${sampleid}.trimmed.fa'
  """
}

process Otu_Preprocess {
  publishDir directoryMap.Otu_preprocess, mode: 'copy', overwrite: true
  tag "$sampleid"
  cpus 32
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
  publishDir directoryMap.Otu_preprocess, mode: 'copy', overwrite: true
  tag "$sampleid"
  cpus 32

  input:
  set val(group), val(sampleid), file(fa), file(count_table), file(align), file(taxonomy) from otu_output

  output:
  set val(group), val(sampleid), file("${sampleid}.trimmed.unique.opti_mcc.list"), file("${sampleid}.trimmed.unique.opti_mcc.0.03.cons.taxonomy"),  file("${sampleid}.otu_table.taxID.filter.dense.dedup.biom"), file("${sampleid}.otu_table.taxID.filter.dense.dedup.txt") into otu_output_similarity97

  script:
  """
  export TERM=xterm
  /app_tools/mothur/mothur "#dist.seqs(processors=${task.cpus}, fasta=${align}, countends=F, cutoff=0.03)" > /dev/null
  /app_tools/mothur/mothur "#cluster(column=${sampleid}.trimmed.unique.dist, count=${count_table}, method=opti, cutoff=0.03)" > /dev/null
  /app_tools/mothur/mothur "#make.shared(list=${sampleid}.trimmed.unique.opti_mcc.list, count=${count_table})" > /dev/null
  /app_tools/mothur/mothur "#classify.otu(list=${sampleid}.trimmed.unique.opti_mcc.list, count=$count_table, cutoff=80, label=0.03, taxonomy=$taxonomy)'" > /dev/null
  awk -F "\t" '{if(substr(\$3,1,7)!="unknown"){print}}' ${sampleid}.trimmed.unique.opti_mcc.0.03.cons.taxonomy > ${sampleid}.trim_unknown.txt

  ## OTUs ref ID, used for comparing with OTU sequence
  /app_tools/mothur/mothur "#make.biom(shared=${sampleid}.trimmed.unique.opti_mcc.shared, constaxonomy=${sampleid}.trimmed.unique.opti_mcc.0.03.cons.taxonomy, matrixtype=dense)" > /dev/null
  mv ${sampleid}.trimmed.unique.opti_mcc.0.03.biom ${sampleid}.otu_table.otuID.biom
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
  """
}

(otu_output_similarity971, otu_output_similarity972) = otu_output_similarity97.into(2)
merge_otu_input = otu_output_similarity971
                .map {
                    it -> it[4]
                }
                .flatten().toList()

process Merge_Otu {
  publishDir directoryMap.Otu_preprocess, mode: 'copy', overwrite: true
  tag "MergedOtu"

  input:
  file(biom) from merge_otu_input

  output:
  set file('otu_table.taxID.filter.dense.dedup.biom'), file('otu_table.taxID.filter.dense.dedup.txt') into merge_otu_output

  script:

  """
  #!/usr/bin/python
  import glob
  import os
  bioms = glob.glob("*.biom")
  list = ','.join(bioms)
  cmd = 'merge_otu_tables.py -i ' + list + ' -o otu_table.taxID.filter.dense.dedup.biom'
  os.system(cmd)
  os.system('. /usr/local/qiime_software/activate.sh')
  cmd = 'biom convert -i otu_table.taxID.filter.dense.dedup.biom -o otu_table.taxID.filter.dense.dedup.txt --to-tsv --header-key taxonomy --table-type="OTU table"'
  os.system(cmd)
  """
}
/*
(merge_otu_output1, merge_otu_output2, merge_otu_output3, merge_otu_output4, merge_otu_output5, merge_otu_output6, merge_otu_output7, merge_otu_output8, merge_otu_output9) = merge_otu_output.into(9)

statistic_analysis_input = otu_output_similarity972   // set val(group), val(sampleid), file(list), file(taxonomy),  file(biom), file(txt)
  .phase(demultiplex_quality_filter_output2) // set val(group), val(sampleid), file("${sampleid}.histograms.txt"), file("${sampleid}.fna")
  .map{
    it->  [it[0][0], it[0][1], it[1][3], it[0][2], it[0][3]]
  }
process Statistic_Analysis {
  publishDir directoryMap.Stat_analysis, mode: 'copy', overwrite: true
  tag "$sampleid"

  input:
  set val(group), val(sampleid), file(fa), file(list),  file(tax) from statistic_analysis_input

  output:
  file("${sampleid}.statistic_sample.txt") into statistic_analysis_output

  script:
  output = sampleid + ".statistic_sample.txt"
  """
  perl /export/EC1680U/perl/bin/16S/statistic_single_sample_mothur.pl ${sampleid} ${fa} ${list} ${tax} ${output}
  """
}
process Statistic_Analysis_part2 {
  publishDir directoryMap.Stat_analysis, mode: 'copy', overwrite: true
  tag "statistic_analysis_part2"

  input:
  file(statistic_sample) from statistic_analysis_output.flatten().toList()
  set file(biom), file(table) from merge_otu_output1

  output:
  set file('Statistic_analysis.png'), file('Statistic_analysis.html'), file('OTU_venn.png') into statistic_analysis_part2_output

  script:

  """
  #!/usr/bin/python
  import glob
  import os
  fns = glob.glob("*.statistic_sample.txt")
  cmd = 'head -1 ' + fns[0] + '> header'
  os.system(cmd)
  cmd = 'rm -rf body'
  os.system(cmd)
  for fn in fns:
    cmd = "sed \'1d\' " + fn + ">> body"
    os.system(cmd)
  cmd = 'cat header body > statistic_sample.txt'
  os.system(cmd)
  cmd = 'plot_StatisticBar.R statistic_sample.txt .'
  os.system(cmd)
  cmd = 'plot_venn.pl ${table} ${mapping_file} OTU_venn.png ${params.group}'
  os.system(cmd)
  """
}
process Alpha_Diversity {
  publishDir directoryMap.alpha_diversity, mode: 'copy', overwrite: true
  tag "alpha_diversity"

  input:
  set file(biom), file(table) from merge_otu_output2

  output:
  set file('otu_table_summary.txt'), file('rarefaction'), file('alpha_div'), file('alpha_div_collated'), file('alpha_rarefaction_plots') into alpha_diversity_output

  script:
  """
  biom summarize-table -i ${biom} -o otu_table_summary.txt
  maximum_reads=\$(grep "Max" otu_table_summary.txt | cut -d" " -f 3 | cut -d"." -f 1)
  if [ \$maximum_reads -ge "30000" ]; then
      multiple_rarefactions.py -i ${biom} -m 1 -x \$maximum_reads -s 2500 -n 10 -o rarefaction/
  else
      multiple_rarefactions.py -i ${biom} -m 1 -x \$maximum_reads -s 250 -n 10 -o rarefaction/
  fi
  alpha_diversity.py -i rarefaction/ -o alpha_div/ -m chao1,enspie,shannon,simpson_reciprocal,fisher_alpha,goods_coverage,observed_species
  collate_alpha.py -i alpha_div/ -o alpha_div_collated/
  make_rarefaction_plots.py -i alpha_div_collated/ -m ${mapping_file} -o alpha_rarefaction_plots
  """
}
process Rank_Abundance {
  publishDir directoryMap.rank_abundance, mode: 'copy', overwrite: true
  tag "alpha_diversity"

  input:
  set file(biom), file(table) from merge_otu_output3

  output:
  set file('rank_abundance_curve.png'), file('rank_abundance_curve_L.png') into rank_abundance_output

  script:
  """
  plot_rank_abundance_graph.py -i ${biom} -s '*' -n -x -f png -o rank_abundance_curve.png
  plot_rank_abundance_graph.py -i ${biom} -s '*'    -x -f png -o rank_abundance_curve_L.png
  """
}
process Alpha_Indices_Table {
  publishDir directoryMap.alpha_indices_table, mode: 'copy', overwrite: true
  tag "alpha_diversity"

  input:
  set file(biom), file(table) from merge_otu_output4

  output:
  file('alpha_indices_table.txt') into alpha_indices_table_output

  script:
  """
  alpha_diversity.py -i ${biom} -t ${params.silva_ref_tree} -m chao1,enspie,shannon,simpson_reciprocal,fisher_alpha,goods_coverage,observed_species -o alpha_indices_table.txt
  """
}
process Taxa_Summary {
  publishDir directoryMap.taxa_summary, mode: 'copy', overwrite: true
  tag "taxa_summary"

  input:
  set file(biom), file(table) from merge_otu_output5

  output:
  set file ('*.biom'), file('*.txt'), file('taxa_summary_plots') into taxa_summary_output

  script:
  """
  echo 'summarize_taxa:level 1,2,3,4,5,6,7' > parameter.txt
  summarize_taxa_through_plots.py -f -o . -i ${biom} -p parameter.txt
  """
}

process Translate {
  publishDir directoryMap.translated, mode: 'copy', overwrite: true
  tag "taxa_summary"

  input:
  set file (biom), file(txt), file('taxa_summary_plots') from taxa_summary_output

  output:
  file ('*') into translate_output

  script:
  """
  translate_Qiime_to_Krnoa.pl otu_table.taxID.filter.dense.dedup_L6.txt .
  ktImportText -o D.1_Krona.html *.txt
  """
}
process Krona {
  publishDir directoryMap.Root_Dir, mode: 'copy', overwrite: true
  tag "taxa_summary"

  input:
  file ('*.txt') from translate_output

  output:
  file ('D.1_Krona.html') into krona_output

  script:
  """
  ktImportText -o D.1_Krona.html *.txt
  """
}
process Heatmap {
  publishDir directoryMap.heatmap, mode: 'copy', overwrite: true
  tag "Heatmap"

  input:
  set file(biom), file(table) from merge_otu_output6

  output:
  file ('*') into heatmap_output

  script:
  """
  plot_heatmap.pl $mapping_file ${table} ${params.group}
  """
}
process Phylogenetic_Tree {
  publishDir directoryMap.phylogenetic_tree, mode: 'copy', overwrite: true
  tag "phylogenetic_tree"

  input:
  set file(biom), file(table) from merge_otu_output7

  output:
  file ('*') into phylogenetic_tree_output

  script:
  """
  /export/EC1680U/perl/bin/16S/plot_phylogenetic_tree.pl ${table} ${params.silva_otu_ref_aligned_fasta}
  """
}

process Metastats {
  publishDir directoryMap.metastats, mode: 'copy', overwrite: true
  tag "Metastats"

  input:
  set file(biom), file(table) from merge_otu_output8

  output:
  file ('*') into metastats_output

  script:
  """
  run_Metastats.pl ${table} $mapping_file . ${params.group}
  """
}

process beta_diversity {
  publishDir directoryMap.beta_diversity, mode: 'copy', overwrite: true
  tag "beta_diversity"

  input:
  set file(biom), file(table) from merge_otu_output9

  output:
  file ('*') into beta_diversity_output

  script:
  """
  echo 'beta_diversity:metrics  euclidean,unweighted_unifrac,weighted_unifrac' > parameter.txt
  beta_diversity_through_plots.py -i ${biom} -f -o . -m $mapping_file -t ${params.silva_ref_tree} -p parameter.txt
  plot_dm_heatmap.R weighted_unifrac_dm.txt weighted_dm_heatmap.png
  plot_dm_heatmap.R unweighted_unifrac_dm.txt unweighted_dm_heatmap.png
  """
}
*/

def defineDirectoryMap(outdir) {
  return [
    'Root_Dir'         : "$outdir/",
    'Qc_RawData'       : "$outdir/B.1_QC/fastqc-rawdata",
    'MQc_RawData'      : "$outdir/B.1_QC/multiqc-rawdata",
    'Trimming_RawData' : "$outdir/B.1_QC/trimming",
    'Qc_TrimedData'    : "$outdir/B.1_QC/fastqc-trimmed",
    'MQc_TrimedData'   : "$outdir/B.1_QC/multiqc-trimmed",
    'Join_Fastq'       : "$outdir/B.2_join_fastq",
    'Dem_Qua_Fil'      : "$outdir/B.3_demultiplex_quilty_filter",
    'Rem_Chimeras'     : "$outdir/B.4_removing_chimeras",
    'Otu_preprocess'   : "$outdir/B.5_OTU_cluster",
    'Stat_analysis'    : "$outdir/B_statistic_analysis",
    'alpha_diversity'  : "$outdir/C_alpha_diversity",
    'rank_abundance'   : "$outdir/C.4_rank_abundance",
    'alpha_indices_table' : "$outdir/C_alpha_indices_table",
    'taxa_summary'     : "$outdir/D.1_taxa_summary",
    'translated'       : "$outdir/D.1_translated",
    'heatmap'          : "$outdir/D.2_heatmap",
    'phylogenetic_tree': "$outdir/D.3_phylogenetic_tree",
    'metastats'        : "$outdir/D.4_metastats",
    'beta_diversity'   : "$outdir/E.1_PCA_PCoA "
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

