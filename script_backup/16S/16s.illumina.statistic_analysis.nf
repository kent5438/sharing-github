#!/usr/bin/env nextflow

/********************************************************************************
 * author: 劉丞原 (daniel.liu)
 * Release date: 2018/08/08
 * Modified date:
 * nextflow run 16s.illumina.nf --rm -c 16.config  2>&1 | tee 16s.illumina.nf.log
 ********************************************************************************/
version = '1.0'

params.file = "../data/Mapping.txt"
params.group = 'Description'
params.lefse_group = 'Description'
if (params.file) {
  mapping_file = file(params.file)
} else {
  mapping_file = file("../data/Mapping.txt")
}

if (params.outdir){
  directoryMap = defineDirectoryMap(params.outdir)
} else{
  exit 1, 'You should specifiy the output directory by using --outdir'
}

// start from Statistic_Analysis

otu_output_similarity97 = Channel
  .from(mapping_file.readLines())
  .map { line ->
    m = line =~ /(.*)\t(.*)\t(.*)\t(.*)\t(.*)/
    id = m[0][1]
    id = id.replaceAll(/-/, "_")
    if (id=~/#SampleID/) {
    } else {
      group = m[0][4]
      list  = "/export/EC1680U/kentchen/16S/MS18119.16S/OUTPUT/B.5_OTU_cluster/97/" + id + ".trimmed.unique.opti_mcc.list"
      tax   = "/export/EC1680U/kentchen/16S/MS18119.16S/OUTPUT/B.5_OTU_cluster/97/" + id + ".trimmed.unique.opti_mcc.0.03.cons.taxonomy"
      [group, id, file(list), file(tax)]
    }
 }

demultiplex_quality_filter_output2 = Channel
  .from(mapping_file.readLines())
  .map { line ->
    m2 = line =~ /(.*)\t(.*)\t(.*)\t(.*)\t(.*)/
    id2 = m2[0][1]
    id2 = id2.replaceAll(/-/, "_")
    if (id2=~/#SampleID/) {
    } else {
      group2 = m2[0][4]
      fa = "/export/EC1680U/kentchen/16S/MS18119.16S/OUTPUT/B.3_demultiplex_quilty_filter/" + id2 + ".fna"
      [group2, id2, file(fa)]
    }
}
(otu_output_similarity971, otu_output_similarity972) = otu_output_similarity97.into(2)

merge_otu_output = Channel.from([file('/export/EC1680U/kentchen/16S/MS18119.16S/OUTPUT/B.5_OTU_cluster/otu_table.taxID.filter.dense.dedup.biom'), file('/export/EC1680U/kentchen/16S/MS18119.16S/OUTPUT/B.5_OTU_cluster/otu_table.taxID.filter.dense.dedup.txt')]).toList()

(merge_otu_output1, merge_otu_output2, merge_otu_output3, merge_otu_output4, merge_otu_output5, merge_otu_output6, merge_otu_output7, merge_otu_output8, merge_otu_output9, merge_otu_output10) = merge_otu_output.into(10)


statistic_analysis_input = otu_output_similarity971
.join(demultiplex_quality_filter_output2, by:[0,1], remainder: true)

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
  //errorStrategy 'ignore' modified by jsyf
  tag "statistic_analysis_part2"

  input:
  file(statistic_sample) from statistic_analysis_output.flatten().toList()
  set file(biom), file(table) from merge_otu_output1

  output:
  set file('Statistic_analysis.png'), file('Statistic_analysis.html'), file('OTU_venn.png') into statistic_analysis_part2_output

  script:

  """
#!/usr/bin/python
import os
fns = []
fh = open("${mapping_file}", "r")
strs = fh.read()
for ln in strs.split('\\n'):
  if not ln.startswith('#'):
    tab = ln.split('\\t')
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
cmd = 'perl /export/EC1680U/perl/bin/16S/plot_venn.pl ${table} ${mapping_file} OTU_venn.png ${params.group}'
os.system(cmd)
  """
}

process Alpha_Diversity {
  publishDir directoryMap.alpha_diversity, mode: 'copy', overwrite: true
  //errorStrategy 'ignore' modified by jsyf
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
  //errorStrategy 'ignore' modified by jsyf
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
  //errorStrategy 'ignore' modified by jsyf
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
process Taxa_Summary {
  publishDir directoryMap.taxa_summary, mode: 'copy', overwrite: true
  //errorStrategy 'ignore' modified by jsyf
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
  perl /export/EC1680U/perl/bin/16S/taxaPlotParser_forCompassPipe.SILVA.pl taxa_summary_plots/raw_data/ taxa_summary_plots/
  """
}
(taxa_summary_output1, taxa_summary_output2, taxa_summary_output3) = taxa_summary_output.into(3)
process Translate {
  publishDir directoryMap.translated, mode: 'copy', overwrite: true
  //errorStrategy 'ignore' modified by jsyf
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
  //errorStrategy 'ignore' modified by jsyf
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
  //errorStrategy 'ignore' modified by jsyf
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
  //errorStrategy 'ignore' modified by jsyf
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
  //errorStrategy 'ignore' modified by jsyf
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
  //errorStrategy 'ignore' modified by jsyf
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
  """
}
process Upgma {
  publishDir directoryMap.upgma, mode: 'copy', overwrite: true
  //errorStrategy 'ignore' modified by jsyf
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
process Lefse {
  publishDir directoryMap.lefse, mode: 'copy', overwrite: true
  //errorStrategy 'ignore' modified by jsyf
  tag "LEfSe"

  input:
  set file (biom), file(txt), file('taxa_summary_plots') from taxa_summary_output2

  output:
  file ('*') into lefsea_output

  script:
  """
  /export/EC1680U/perl/bin/16S/run_LEfSe.pl otu_table.taxID.filter.dense.dedup_L6.txt ${mapping_file} ${params.lefse_group} .
  """
}
(beta_diversity_output1, beta_diversity_output2,beta_diversity_output3) = beta_diversity_output.into(3)
process Nmds {
  publishDir directoryMap.nmds, mode: 'copy', overwrite: true
  //errorStrategy 'ignore' modified by jsyf
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
  """
}
process Anosim {
  publishDir directoryMap.anosim, mode: 'copy', overwrite: true
  //errorStrategy 'ignore' modified by jsyf
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
  //errorStrategy 'ignore' modified by jsyf
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
  //errorStrategy 'ignore' modified by jsyf
  tag "DCA/RDA/CCA"

  input:
  set file (biom), file(txt), file('taxa_summary_plots') from taxa_summary_output3

  output:
  file ('*') into dca_rda_cca_output

  script:
  """
  sed 's/#/_/g' otu_table.taxID.filter.dense.dedup_L7.txt > otu_table.taxID.filter.dense.dedup_L7.txt.1
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
