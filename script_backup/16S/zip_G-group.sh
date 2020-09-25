#! /bin/bash

for i in OUTPUT_Merge_G*; do cd ${i}/; /export/EC1680U/perl/bin/16S/g16S_report.sh Illu; cd ..; done

for i in `ls -d OUTPUT_Merge_G*/ | cut -d'/' -f1 | cut -d'_' -f3`; do mv OUTPUT_Merge_${i}/report/_reports OUTPUT_Merge_${i}/report/MS19221-${i}_NPUST-JHShyu_report_191209; zip -r OUTPUT_Merge_${i}/report/MS19221-${i}_NPUST-JHShyu_report_191209 OUTPUT_Merge_${i}/report/MS19221-${i}_NPUST-JHShyu_report_191209; done
