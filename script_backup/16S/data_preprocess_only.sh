#!/bin/bash

# basic parameter
ID=`id -u -n`
workPath=$PWD
mapping_file="$workPath/data/Mapping.txt" # absolute path better, files must be fastq format.
group_title="Description"
venn_group="$group_title"
lefse_group="$group_title"

# FLASH v1.2.11 parameter
min_overlap=10       # -m,  The minimum required overlap length between two reads to provide a confident overlap.  Default: 10bp.
max_overlap=65       # -M,  Maximum overlap length expected in approximately 90% of read pairs.  It is by default set to 65bp, which works well for 100bp reads generated from a 180bp library.
mismatch_ratio=0.25  # -x,  Maximum allowed ratio between the number of mismatched base pairs and the overlap length. Default: 0.25.

# AlienTrimmer v0.4 parameter
trim_fasta="$workPath/data/primer.fa" # -c,  for single-end, pair-end (absolute path better, leave '-' for none)
trim_forward_fasta="-"                # -cf, for pair-end only (absolute path better, leave '-' for none)
trim_reverse_fasta="-"                # -cr, for pair-end only (absolute path better, leave '-' for none)
kmer=10;         # -k, k-mer decomposition, between 5 and 15, default 10
mismatches=5;    # -m, maximum mismatches, between 0 and 15, default kmer/2
read_length=150; # -l, minimum read length, default 15
phred_quality=1; # -q, Phred quality score cut-off, between 0 and 40, default 20
percentage=50;   # -p, minimum allowed percentage of correctly called nucleotides, between 0 and 100, default 0

# Mothur parameters
classify_cutoff=60

# other parameters
phred_offset=33
quality=19
num_permutations=999
sequence_depth=110
thread=32

# reference
chimeras_ref='/references/gold.fa'
otu_ref_path='/references/gg_13_8_otus'
otu_ref_tree="$otu_ref_path/trees/97_otus.tree"
mothur_otu_ref_path='/references/gg_13_8_99'
mothur_otu_ref_map="$otu_ref_path/otus/99_otu_map.txt"
mothur_otu_ref_fasta="$mothur_otu_ref_path/gg_13_8_99.fasta"
mothur_otu_ref_aligned_fasta="$mothur_otu_ref_path/gg_13_8_99.refalign"
mothur_otu_ref_taxonomy="$mothur_otu_ref_path/gg_13_8_99.gg.tax"

# env setting
. /usr/local/qiime_software/activate.sh
export PATH=$PATH:/opt/local/software/fastx_toolkit_0.0.13/bin
export AlienTrimmer='/app_tools/AlienTrimmer.jar'
export PATH=/app_tools/mothur/:$PATH



#-- Start --#
START=$(date +%s)

echo 'B data preprocess'
$workPath/data_preprocess.pl $mapping_file $phred_offset $quality $min_overlap $max_overlap $mismatch_ratio $trim_fasta $trim_forward_fasta $trim_reverse_fasta $kmer $mismatches $read_length $phred_quality $percentage $thread '    '

