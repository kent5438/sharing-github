#!/usr/bin/env nextflow

/********************************************************************************
 * author: 劉丞原 (daniel.liu)
 * Release date: 2018/08/23
 * Modified date:
 * nextflow run 16s.illumina.nf --rm -c 16.config  2>&1 | tee 16s.illumina.nf.log
 ********************************************************************************/
version = '1.0'

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
    m = line =~ /(.*)\t(.*)\t(.*)\t(.*)\t(.*)/
    id = m[0][1]
    id = id.replaceAll(/-/, "_")
    if (id=~/#SampleID/) {
    } else {
      desc = m[0][4]
      reads = m[0][5]
      [desc, id, reads]
    }
   }

(sample_reads_qc, sample_reads_trim) = sample_reads.into(2)
process Qc_Raw_Data {
  publishDir directoryMap.qc_rawData, mode: 'copy', overwrite: true
  tag "$sampleid"

  input:
  set val(group), val(sampleid), val(reads) from sample_reads_qc

  output:
  set val(group), val(sampleid), file('*_fastqc.{zip,html}') into fastqc_rawdata_output

  script:
  """
  fastqc -o . ${reads}
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
  multiqc -f -q -c /export/EC1680U/perl/bin/16S/PB-16S_multiqc_config.yaml -o . .
  """
}
process Trimming_Raw_Data {
  publishDir directoryMap.trimming_rawdata, mode: 'copy', overwrite: true
  tag "$sampleid"

  input:
  set val(group), val(sampleid), val(reads) from sample_reads_trim

  output:
  set val(group), val(sampleid), file(trim) into trim_rawdata_output

  script:
  trim = sampleid + ".clean.fastq"

  """
  cutadapt -m 600 --discard-untrimmed -g AGRGTTYGATYMTGGCTCAG...AAGTCGTAACAAGGTARCY -g RGYTACCTTGTTACGACTT...CTGAGCCAKRATCRAACYCT -o ${trim} ${reads} > /dev/null
  """
}
(trim_rawdata_output_qc, trim_rawdata_output_dem_qua_filter) = trim_rawdata_output.into(2)
process Qc_Trimed_Data {
  publishDir directoryMap.qc_trimeddata, mode: 'copy', overwrite: true
  tag "$sampleid"

  input:
  set val(group), val(sampleid), file(trim) from trim_rawdata_output_qc

  output:
  set val(group), val(sampleid), file('*_fastqc.{zip,html}') into fastqc_trimdedata_output

  script:
  """
  fastqc -o . ${trim}
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
  multiqc -f -q -c /export/EC1680U/perl/bin/16S/PB-16S_multiqc_config.yaml -o . .
  """
}
process Dem_Qua_Filter {
  publishDir directoryMap.dem_qua_fil, mode: 'copy', overwrite: true
  tag "$sampleid"

  input:
  set val(group), val(sampleid), file(fastq) from trim_rawdata_output_dem_qua_filter

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
#  /opt/mothur_v1.44.0/mothur "#chimera.uchime(processors=${task.cpus}, fasta=${fna}, reference=${params.chimeras_ref})" > /dev/null 2>&1
  /export/EC1680U/kentchen/opt/usearch11.0.667_i86linux32 -uchime2_ref ${fna} -db ${params.chimeras_ref} -strand plus -threads ${task.cpus} -mode high_confidence -chimeras ${sampleid}.ref.uchime.chimeras -notmatched ${sampleid}.trimmed.fa
#  removing_chimera.pl '${fna}' '${sampleid}.ref.uchime.chimeras' '${sampleid}.trimmed.fa'
  """
}

(removing_chimeras_output1,removing_chimeras_output2) = removing_chimeras_output.into(2)

process Otu_Preprocess {
  publishDir directoryMap.otu_preprocess, mode: 'copy', overwrite: true
  tag "$sampleid"
  cpus 20

  input:
  set val(group), val(sampleid), file(fa) from removing_chimeras_output1

  output:
  set val(group), val(sampleid), file("${sampleid}.trimmed.unique.fa"), file("${sampleid}.trimmed.count_table"), file("${sampleid}.trimmed.unique.align"),  file("${sampleid}.trimmed.unique.silva.wang.taxonomy") into otu_output

  script:
  """
  export TERM=xterm
  mothur_make_group.pl ${fa}
  /opt/mothur_v1.44.0/mothur "#unique.seqs(fasta=${fa})" > /dev/null # 輸出 XX.trimmed.unique.fa & XX.trimmed.names
  /opt/mothur_v1.44.0/mothur "#count.seqs(name=${sampleid}.trimmed.names, group=group.txt)" > /dev/null # 輸出 XX.trimmed.count_table
  /opt/mothur_v1.44.0/mothur "#classify.seqs(processors=${task.cpus}, fasta=${sampleid}.trimmed.unique.fa, name=${sampleid}.trimmed.names, template=${params.silva_otu_ref_fasta}, taxonomy=${params.silva_otu_ref_taxonomy}, cutoff=80, method=wang)"  # 輸出 XX.trimmed.unique.raw_taxonomy.wang.taxonomy  & XX.trimmed.unique.raw_taxonomy.wang.tax.summary
  /opt/mothur_v1.44.0/mothur "#remove.lineage(fasta=${sampleid}.trimmed.unique.fa, count=${sampleid}.trimmed.count_table, taxonomy=${sampleid}.trimmed.unique.raw_taxonomy.wang.taxonomy, taxon='o__Chloroplast-f__Mitochondria-unknown')" > /dev/null
  mv ${sampleid}.trimmed.unique.raw_taxonomy.wang.pick.taxonomy ${sampleid}.trimmed.unique.silva.wang.taxonomy
  if [ -f "${sampleid}.trimmed.unique.pick.fa" ]; then
    mv ${sampleid}.trimmed.unique.pick.fa ${sampleid}.trimmed.unique.fa
  fi
  if [ -f "${sampleid}.trimmed.unique.pick.count_table" ]; then
    mv ${sampleid}.trimmed.pick.count_table ${sampleid}.trimmed.count_table
  fi
  /opt/mothur_v1.44.0/mothur "#align.seqs(processors=${task.cpus}, candidate=${sampleid}.trimmed.unique.fa, template=${params.silva_otu_ref_aligned_fasta})" > /dev/null # 輸出 XX.trimmed.unique.align & XX.trimmed.unique.align.report && XX.trimmed.unique.flip.accnos
  """
}

process Otu_Similarity97 {
  publishDir directoryMap.otu_preprocess, mode: 'copy', overwrite: true
  tag "$sampleid"
  cpus 20

  input:
  set val(group), val(sampleid), file(fa), file(count_table), file(align), file(taxonomy) from otu_output

  output:
  set val(group), val(sampleid), file("97/${sampleid}.trimmed.unique.opti_mcc.list"), file("97/${sampleid}.trimmed.unique.opti_mcc.0.01.cons.taxonomy"),  file("97/${sampleid}.otu_table.taxID.filter.dense.dedup.biom"), file("97/${sampleid}.otu_table.taxID.filter.dense.dedup.txt"), file('97/*') into otu_output_similarity97

  script:
  """
  export TERM=xterm
  mkdir 97
  cp ${fa} 97
  cp ${count_table} 97
  cp ${align} 97
  cp ${taxonomy} 97
  cd 97/


  /opt/mothur_v1.44.0/mothur "#dist.seqs(processors=${task.cpus}, fasta=${align}, countends=F, cutoff=0.01)" > /dev/null
  /opt/mothur_v1.44.0/mothur "#cluster(column=${sampleid}.trimmed.unique.dist, count=${count_table}, method=opti, cutoff=0.01)" > /dev/null
  /opt/mothur_v1.44.0/mothur "#make.shared(list=${sampleid}.trimmed.unique.opti_mcc.list, count=${count_table})" > /dev/null
  /opt/mothur_v1.44.0/mothur "#classify.otu(list=${sampleid}.trimmed.unique.opti_mcc.list, count=${count_table}, cutoff=80, label=0.01, taxonomy=$taxonomy)'" > /dev/null

  /opt/mothur_v1.44.0/mothur "#remove.lineage(constaxonomy=${sampleid}.trimmed.unique.opti_mcc.0.01.cons.taxonomy, shared=${sampleid}.trimmed.unique.opti_mcc.shared, taxon='unknown', label=0.01)"
  mv ${sampleid}.trimmed.unique.opti_mcc.0.01.cons.pick.cons.taxonomy ${sampleid}.trimmed.unique.opti_mcc.0.01.cons.taxonomy
  if [ -e "${sampleid}.trimmed.unique.opti_mcc.0.01.pick.shared" ]; then
    mv ${sampleid}.trimmed.unique.opti_mcc.0.01.pick.shared ${sampleid}.trimmed.unique.opti_mcc.shared
  fi
  
## OTUs ref ID, used for comparing with OTU sequence
  /opt/mothur_v1.44.0/mothur "#make.biom(shared=${sampleid}.trimmed.unique.opti_mcc.shared, constaxonomy=${sampleid}.trimmed.unique.opti_mcc.0.01.cons.taxonomy, matrixtype=dense)" > /dev/null
  mv ${sampleid}.trimmed.unique.opti_mcc.0.01.biom ${sampleid}.otu_table.otuID.biom


  filter_otus_from_otu_table.py -i ${sampleid}.otu_table.otuID.biom -o ${sampleid}.otu_table.otuID.filter.biom --min_count_fraction 0.000000000000001 # ref: https://goo.gl/EMwR9N
  biom convert -i ${sampleid}.otu_table.otuID.filter.biom -o ${sampleid}.otu_table.otuID.filter.txt --to-tsv --header-key taxonomy --table-type="OTU table"

  ## Generate OTUs sequence and taxonomy output
  /opt/mothur_v1.44.0/mothur "#get.oturep(column= ${sampleid}.trimmed.unique.dist, list= ${sampleid}.trimmed.unique.opti_mcc.list, count=${count_table}, fasta=${fa}, method=abundance, cutoff=0.01)" > /dev/null

  ## Convert otuID.filter to taxID.filter
#  vsearch --usearch_global  ${sampleid}.trimmed.unique.opti_mcc.0.01.rep.fasta --db ${params.silva_otu_ref_fasta} --id 0.1 --iddef 1 --blast6out results.blast6 --maxaccepts 1 --threads 32
#  blastn -query ${sampleid}.trimmed.unique.opti_mcc.0.01.rep.fasta -db ${params.silva_otu_ref_fasta} -num_threads 32 -outfmt 6 -max_target_seqs 1 -max_hsps 1 -out ${sampleid}.results.blast6

  ### Contributed by JSYF ###
  perl /export/EC1680U/perl/bin/16S/biomProc.pl sparseToDense \
  	${sampleid}.otu_table.otuID.filter.biom \
  	${sampleid}.otu_table.otuID.filter.dense.biom

#  perl /export/EC1680U/perl/bin/16S/biomProc.pl idChangeUsingFasta \
#  	${sampleid}.otu_table.otuID.filter.dense.biom \
#  	${sampleid}.otu_table.taxID.filter.dense.biom \
#  	${sampleid}.trimmed.unique.opti_mcc.0.01.rep.fasta,${sampleid}.results.blast6
  python /export/EC1680U/perl/bin/16S/otuToTax.py \
    ${sampleid}.otu_table.otuID.filter.dense.biom \
    ${sampleid}.otu_table.taxID.filter.dense.biom

  perl /export/EC1680U/perl/bin/16S/biomProc.pl uniqById \
  	${sampleid}.otu_table.taxID.filter.dense.biom \
  	${sampleid}.otu_table.taxID.filter.dense.dedup.biom
  ########################

  biom convert -i ${sampleid}.otu_table.taxID.filter.dense.dedup.biom -o ${sampleid}.otu_table.taxID.filter.dense.dedup.txt --to-tsv --header-key taxonomy --table-type="OTU table"
  rm ${fa}
  rm ${count_table}
  rm ${align}
  rm ${taxonomy}
#  rm results.blast6
  rm mothur.*.logfile
  cd ../
  """
}

(otu_output_similarity971, otu_output_similarity972) = otu_output_similarity97.into(2)
merge_otu_input = otu_output_similarity971
                .map {
                    it -> it[4]
                }
                .flatten().toList()

process Merge_Otu {
  publishDir directoryMap.otu_preprocess, mode: 'copy', overwrite: true
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
bioms = sorted(glob.glob("*.biom"))
list = ','.join(bioms)
cmd = "cut -f1 ${mapping_file} | grep -v '#' > sample_id.txt"
os.system(cmd)
cmd = 'merge_otu_tables.py -i ' + list + ' -o otu_table.taxID.filter.dense.dedup.biom'
os.system(cmd)
cmd = "sort_otu_table.py -i otu_table.taxID.filter.dense.dedup.biom -o tmp -l sample_id.txt"
os.system(cmd)
cmd = "mv tmp otu_table.taxID.filter.dense.dedup.biom"
os.system(cmd)
cmd = 'biom convert -i otu_table.taxID.filter.dense.dedup.biom -o otu_table.taxID.filter.dense.dedup.txt --to-tsv --header-key taxonomy --table-type="OTU table"'
os.system(cmd)
  """
}
(merge_otu_output1, merge_otu_output2, merge_otu_output3, merge_otu_output4, merge_otu_output5, merge_otu_output6, merge_otu_output7, merge_otu_output8, merge_otu_output9, merge_otu_output10, merge_otu_output11) = merge_otu_output.into(11)

// set val(group), val(sampleid), file(list), file(taxonomy),  file(biom), file(txt), file('97')
tmp_otu_output_similarity972 = otu_output_similarity972
.map {
  it->  [it[0], it[1], it[2], it[3]]
}
.collect()
.flatten()
.collate(4)
//.view()

// set val(group), val(sampleid), file("${sampleid}.fna")
tmp_removing_chimeras_output2 = removing_chimeras_output2
.collect()
.flatten()
.collate(3)
//.view()

statistic_analysis_input = tmp_otu_output_similarity972
.join(tmp_removing_chimeras_output2, by:[0,1], remainder: true)
//.view()

process Statistic_Analysis {
  publishDir directoryMap.stat_analysis, mode: 'copy', overwrite: true
  tag "$sampleid"

  input:
  set val(group), val(sampleid), file(list), file(tax), file(fa) from statistic_analysis_input

  output:
  file("${sampleid}.statistic_sample.txt") into statistic_analysis_output

  script:
  output = sampleid + ".statistic_sample.txt"
  """
  perl /export/EC1680U/perl/bin/16S/statistic_single_sample_mothur.pl ${sampleid} ${fa} ${list} ${tax} ${output}
  """
}

process Statistic_Analysis_part2 {
  publishDir directoryMap.stat_analysis, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "statistic_analysis_part2"

  input:
  file(statistic_sample) from statistic_analysis_output.flatten().toList()

  output:
  set file('Statistic_analysis.png'), file('Statistic_analysis.html') into statistic_analysis_part2_output

  script:
  """
  #!/usr/bin/python
  import os
  fns = []
  fh = open('${mapping_file}', "r")
  strs = fh.read()
  for ln in strs.split('\\n'):
    if not ln.startswith('#'):
      tab = ln.split('\t')
      fn = tab[0] + ".statistic_sample.txt"
      fns.append(fn)
  fh.close()

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
  """
}

// 逸帆修改版本
process Alpha_Diversity {
  publishDir directoryMap.alpha_diversity, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
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

  LINE=`grep -n 'If the lines' alpha_rarefaction_plots/rarefaction_plots.html | awk -F':' '{print \$1}'`
  NEXT_LINE=\$((\$LINE+1))
  STRING='<br><br><span style="color: black; font-family:Arial,Verdana; font-size:14; font-weight:bold;"><a href="./rarePlot_data.xls" target="_blank"> View full table (.xls) </a></span>'
  sed -i "\${NEXT_LINE}i \$STRING" alpha_rarefaction_plots/rarefaction_plots.html
  perl /export/EC1680U/perl/bin/16S/rarePlotParser_forCompassPipe.pl alpha_div_collated/ alpha_rarefaction_plots/
  """
}
process Rank_Abundance {
  publishDir directoryMap.rank_abundance, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "rank_abundance"

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
  errorStrategy 'ignore'
  tag "alpha_indices_table"

  input:
  set file(biom), file(table) from merge_otu_output4

  output:
  file('alpha_indices_table.txt') into alpha_indices_table_output

  script:
  """
  alpha_diversity.py -i ${biom} -t ${params.silva_ref_tree} -m chao1,enspie,shannon,simpson_reciprocal,fisher_alpha,goods_coverage,observed_species -o alpha_indices_table.txt
  """
}
// 逸帆修改版本
process Taxa_Summary {
  publishDir directoryMap.taxa_summary, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "new_taxa_summary"

  input:
  set file(biom), file(table) from merge_otu_output5

  output:
  set file ('*.biom'), file('*.txt'), file('taxa_summary_plots') into taxa_summary_output

  script:

  """
  summarize_taxa.py -i ${biom} -o taxa_summary_plots -a -L 1,2,3,4,5,6,7

  plot_taxa_summary.py \
          -i taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L1.txt,taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L2.txt,taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L3.txt,taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L4.txt,taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L5.txt,taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L6.txt,taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L7.txt \
          -l Kingdom,Phylum,Class,Order,Family,Genus,Species \
          -c area,bar -o taxa_summary_plots

  cd taxa_summary_plots
  for i in otu_table.taxID.filter.dense.dedup_L*.txt
  do
          sed -i '1d' \$i
          sed -i 's/^#OTU ID/Taxon/' \$i
  done
  cd ..

  sed -i '/Kingdom/i<tr><td class="normal" colspan=2 style="font-size:14"><a href="./taxaPlot_data.xls" target="_blank"> View full table (.xls) </a></td></tr>' taxa_summary_plots/area_charts.html
  sed -i '/Kingdom/i<tr><td class="normal" colspan=2 style="font-size:14"><a href="./taxaPlot_data.xls" target="_blank"> View full table (.xls) </a></td></tr>' taxa_summary_plots/bar_charts.html
  mv -f taxa_summary_plots/otu_table.taxID.filter.dense.dedup_L* taxa_summary_plots/raw_data/
  cp taxa_summary_plots/raw_data/* .
  perl /export/EC1680U/perl/bin/16S/taxaPlotParser_forCompassPipe.SILVA.xlsx.pl taxa_summary_plots/raw_data/ taxa_summary_plots/
  """
}

(taxa_summary_output1, taxa_summary_output2, taxa_summary_output3) = taxa_summary_output.into(3)
process Translate {
  publishDir directoryMap.translated, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "translate"

  input:
  set file (biom), file(txt), file('taxa_summary_plots') from taxa_summary_output1

  output:
  file ('*') into translate_output

  script:
  """
  /export/EC1680U/perl/bin/16S/translate_Qiime_to_Krnoa.pl otu_table.taxID.filter.dense.dedup_L6.txt .
  """
}
process Krona {
  publishDir directoryMap.root_dir, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "krona"

  input:
  file (txt) from translate_output

  output:
  file ('D.1_Krona.html') into krona_output

  script:
  """
  ktImportText -o D.1_Krona.html *.txt
  """
}
process Heatmap {
  publishDir directoryMap.heatmap, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "heatmap"

  input:
  set file(biom), file(table) from merge_otu_output6

  output:
  file ('*') into heatmap_output

  script:
  """
  plot_heatmap.pl ${mapping_file} ${table} ${params.group}
  """
}
process Phylogenetic_Tree {
  publishDir directoryMap.phylogenetic_tree, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "phylogenetic_tree"

  input:
  set file(biom), file(table) from merge_otu_output7

  output:
  file ('*') into phylogenetic_tree_output

  script:
  """
  export TERM=xterm
  /export/EC1680U/perl/bin/16S/plot_phylogenetic_tree.pl ${table} ${params.silva_otu_ref_aligned_fasta}
  """
}
process Metastats {
  publishDir directoryMap.metastats, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "metastats"

  input:
  set file(biom), file(table) from merge_otu_output8

  output:
  file ('*') into metastats_output

  script:
  """
  #!/usr/bin/python
  import glob
  import os

  def unique(list1):
    unique_list = []
    for x in list1:
      if x not in unique_list:
        unique_list.append(x)
    return unique_list

  with open('${mapping_file}', 'r') as fn:
    lines = fn.read().splitlines()

  columns = lines[0].split('\t')
  index = [i for i, elem in enumerate(columns) if '${params.group}' in elem]
  groups = []
  for idx in range(1, len(lines)):
    line = lines[idx]
    columns = line.split('\t')
    groups.append(columns[index[0]])

  UG = unique(groups)
  if (len(UG) == 2):
    cmd = 'run_Metastats.pl ${table} ${mapping_file} . ${params.group}'
    os.system(cmd)
  """
}

process Beta_Diversity {
  publishDir directoryMap.beta_diversity, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "beta_diversity"

  input:
  set file(biom), file(table) from merge_otu_output9

  output:
  file ('*') into beta_diversity_output

  script:
  """
  echo 'beta_diversity:metrics  euclidean,unweighted_unifrac,weighted_unifrac' > parameter.txt
  beta_diversity_through_plots.py -i ${biom} -f -o . -m ${mapping_file} -t ${params.silva_ref_tree} -p parameter.txt
  plot_dm_heatmap.R weighted_unifrac_dm.txt weighted_dm_heatmap.png
  plot_dm_heatmap.R unweighted_unifrac_dm.txt unweighted_dm_heatmap.png
# Plot 2D PCoA (add by kentchen)
#  make_2d_plots.py -i weighted_unifrac_pc.txt -m ${mapping_file} -o weighted_unifrac_2d_pcoa_plots
  /export/EC1680U/kentchen/miniconda3/envs/phylotoast/bin/PCoA.py -i weighted_unifrac_pc.txt -m ${mapping_file} -b Description -o weighted_unifrac_pc_2D_PCoA_plots.png --annotate_points -t "Weighted PCoA 2D plot" --figsize 12 6
  """
}
process Upgma {
  publishDir directoryMap.upgma, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "upgma"
  cpus 20
  input:
  set file(biom), file(table) from merge_otu_output10

  output:
  set file ('rarefaction'), file('weighted_unifrac'), file('unweighted_unifrac') into upgma_output

  script:
  """
  #  E.2 Clustering analysis of UPGMA'
  mkdir rarefaction
  mkdir -p weighted_unifrac/rare_dm
  multiple_rarefactions_even_depth.py -i otu_table.taxID.filter.dense.dedup.biom -d 110 -o rarefaction
  #  weighted
  parallel_beta_diversity.py -i rarefaction/ -o weighted_unifrac/rare_dm/ -m weighted_unifrac -t /references/SILVA_132_release/trees/99/99_otus.tre -O 20
  rm -rf weighted_unifrac/rare_dm/BDIV*
  mkdir -p weighted_unifrac/rare_upgma
  upgma_cluster.py -i weighted_unifrac/rare_dm/ -o weighted_unifrac/rare_upgma/
  consensus_tree.py -i weighted_unifrac/rare_upgma/ -o weighted_unifrac/rare_upgma_consensus.tre
  mkdir -p weighted_unifrac/upgma_cmp
  tree_compare.py -s weighted_unifrac/rare_upgma/ -m weighted_unifrac/rare_upgma_consensus.tre -o weighted_unifrac/upgma_cmp/
  /export/EC1680U/perl/bin/16S/plot_tree.pl weighted_unifrac/upgma_cmp/master_tree.tre ${mapping_file} weighted_unifrac/weighted_UPGMA_tree.png Description "Weighted UPGMA tree"
  # unweighted
  mkdir -p unweighted_unifrac/rare_dm
  parallel_beta_diversity.py -i rarefaction/ -o unweighted_unifrac/rare_dm/ -m unweighted_unifrac -t /references/SILVA_132_release/trees/99/99_otus.tre -O 20
  rm -rf unweighted_unifrac/rare_dm/BDIV*
  mkdir -p unweighted_unifrac/rare_upgma
  upgma_cluster.py -i unweighted_unifrac/rare_dm -o unweighted_unifrac/rare_upgma/
  consensus_tree.py -i unweighted_unifrac/rare_upgma -o unweighted_unifrac/rare_upgma_consensus.tre
  mkdir -p unweighted_unifrac/upgma_cmp
  tree_compare.py -s unweighted_unifrac/rare_upgma -m unweighted_unifrac/rare_upgma_consensus.tre -o unweighted_unifrac/upgma_cmp/
  /export/EC1680U/perl/bin/16S/plot_tree.pl unweighted_unifrac/upgma_cmp/master_tree.tre ${mapping_file} unweighted_unifrac/unweighted_UPGMA_tree.png Description "Unweighted UPGMA tree"
  """
}


process Tax4Fun {
  publishDir directoryMap.tax4fun, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "tax4fun"
  container "daniel0325/tax4fun:0.0.0.1"

  input:
  set file(biom), file(table) from merge_otu_output11

  output:
  set file ('Kegg_pathway2*'), file('tax4fun*'), file('pathway.stat*') into tax4fun_output

  script:
  """
  cat $table | sed 's/[a-z]__//g' > motu_table.taxID.filter.dense.dedup.txt
  awk -F "\t" '{print \$1"\t"\$4}' ${mapping_file} > group.txt
  /opt/conda/bin/Rscript /export/EC1680U/perl/bin/16S/Tax4Fun/Tax4Fun.r
  """
}

process Lefse {
  publishDir directoryMap.lefse, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "LEfSe"

  input:
  set file (biom), file(txt), file('taxa_summary_plots') from taxa_summary_output2

  output:
  file ('*') into lefsea_output

  script:
  """
#  /export/EC1680U/perl/bin/16S/run_LEfSe.pl otu_table.taxID.filter.dense.dedup_L6.txt ${mapping_file} ${params.lefse_group} .
  perl /export/EC1680U/perl/bin/16S/LEfSe_preprocessInNextflow.pl otu_table.taxID.filter.dense.dedup_L6.txt
  /export/EC1680U/perl/bin/16S/run_LEfSe.pl otu_table.taxID.filter.dense.dedup_L6.charact.txt ${mapping_file} ${params.lefse_group} .
  """
}
(beta_diversity_output1, beta_diversity_output2,beta_diversity_output3) = beta_diversity_output.into(3)
process Nmds {
  publishDir directoryMap.nmds, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "NMDS"

  input:
  file ('*') from beta_diversity_output1

  output:
  file('*') into nmds_output

  script:
  """
  mkdir -p PCA
  mkdir -p PCoA_weighted
  mkdir -p PCoA_unweighted
  nmds.py -i euclidean_dm.txt -o PCA_nmds.txt
  convert_nmds_to_emperor.py PCA_nmds.txt PCA_emperor.txt
  make_emperor.py -i PCA_emperor.txt --map_fp ${mapping_file} -o PCA

  nmds.py -i weighted_unifrac_dm.txt -o PCoA_weighted_nmds.txt
  convert_nmds_to_emperor.py PCoA_weighted_nmds.txt PCoA_weighted_emperor.txt
  make_emperor.py -i PCoA_weighted_emperor.txt --map_fp ${mapping_file} -o PCoA_weighted

  nmds.py -i unweighted_unifrac_dm.txt -o PCoA_unweighted_nmds.txt
  convert_nmds_to_emperor.py PCoA_unweighted_nmds.txt PCoA_unweighted_emperor.txt
  make_emperor.py -i PCoA_unweighted_emperor.txt --map_fp ${mapping_file} -o PCoA_unweighted
# Plot 2D NMDS (add by kentchen)
  perl /export/EC1680U/perl/bin/16S/nmdsMerge.pl ${mapping_file} PCoA_weighted_nmds.txt
  Rscript /export/EC1680U/Rscript/16S/nmds_plot.R
  """
}
process Anosim {
  publishDir directoryMap.anosim, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "Anosim"

  input:
  file ('*') from beta_diversity_output2

  output:
  set file ('PCA'), file('PCoA_weighted'), file('PCoA_unweighted') into anosim_output

  script:
  """
  mkdir -p PCA
  mkdir -p PCoA_weighted
  mkdir -p PCoA_unweighted
  if [ \$(group_status.pl ${mapping_file} "ANOSIM" ${params.group}) -eq 1 ]
  then
      compare_categories.py -i euclidean_dm.txt -o PCA -m ${mapping_file} -c ${params.group} --method anosim --num_permutations ${params.num_permutations}
      compare_categories.py -i weighted_unifrac_dm.txt -o PCoA_weighted -m ${mapping_file} -c ${params.group} --method anosim --num_permutations ${params.num_permutations}
      compare_categories.py -i unweighted_unifrac_dm.txt -o PCoA_unweighted -m ${mapping_file} -c ${params.group} --method anosim --num_permutations ${params.num_permutations}
  fi
  """
}
process Mrpp {
  publishDir directoryMap.mrpp, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "MRPP"

  input:
  file ('*') from beta_diversity_output3

  output:
  set file ('PCA'), file('PCoA_weighted'), file('PCoA_unweighted') into mrpp_output

  script:
  """
  mkdir -p PCA
  mkdir -p PCoA_weighted
  mkdir -p PCoA_unweighted
  if [ \$(group_status.pl ${mapping_file} "MRPP" "${params.group}") -eq 1 ]
  then
      compare_categories.py -i euclidean_dm.txt -o PCA -m ${mapping_file} -c ${params.group} --method mrpp --num_permutations ${params.num_permutations}
      compare_categories.py -i weighted_unifrac_dm.txt -o PCoA_weighted -m ${mapping_file} -c ${params.group} --method mrpp --num_permutations ${params.num_permutations}
      compare_categories.py -i unweighted_unifrac_dm.txt -o PCoA_unweighted -m ${mapping_file} -c ${params.group} --method mrpp --num_permutations ${params.num_permutations}
  fi
  """
}
process dca_rda_cca {
  publishDir directoryMap.dca_rda_cca, mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  tag "DCA/RDA/CCA"

  input:
  set file (biom), file(txt), file('taxa_summary_plots') from taxa_summary_output3

  output:
  file ('*') into dca_rda_cca_output

  script:
  """
  sed 's/#/_/g' otu_table.taxID.filter.dense.dedup_L7.txt > otu_table.taxID.filter.dense.dedup_L7.txt.1
#  plot_3A.pl otu_table.taxID.filter.dense.dedup_L7.txt.1 ${mapping_file} 'Species' .
  /export/EC1680U/perl/bin/16S/plot_3A.pl otu_table.taxID.filter.dense.dedup_L7.txt.1 ${mapping_file} 'Species' .
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
