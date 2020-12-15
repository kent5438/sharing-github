#! /usr/bin/env python

### % . activate qiime ###
### % locate at 97     ###

import glob
import os
bioms = glob.glob("*.otu_table.taxID.filter.dense.dedup.biom")
list = ','.join(bioms)
#cmd = "cut -f1 /export/EC1680U/kentchen/16S/MS18182-MS19042-MS19182.16S/mouse-re/data/Mapping.txt | grep -v '#' > sample_id.txt"
cmd = "cut -f1 ../../../data/Mapping_all.txt | grep -v '#' > sample_id.txt"
os.system(cmd)
cmd = 'merge_otu_tables.py -i ' + list + ' -o otu_table.taxID.filter.dense.dedup.biom'
os.system(cmd)
cmd = "sort_otu_table.py -i otu_table.taxID.filter.dense.dedup.biom -o tmp -l sample_id.txt"
os.system(cmd)
cmd = "mv tmp otu_table.taxID.filter.dense.dedup.biom"
os.system(cmd)
cmd = 'biom convert -i otu_table.taxID.filter.dense.dedup.biom -o otu_table.taxID.filter.dense.dedup.txt --to-tsv --header-key taxonomy --table-type="OTU table"'
os.system(cmd)
