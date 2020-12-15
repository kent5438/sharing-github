#! /bin/bash
# Author: kentchen
# Date: 2017-06-12

DIR=`pwd`
echo "* You're now at: $DIR"

VAR=0

while [ $# -gt 0 ]
do
OPT=$1

case $OPT in
    -h|--help)
    HELP=$1
    ((VAR++))
    ;;
    -i|--input)
    INPUT=$2
    shift
    ((VAR++))
    ;;
    -c|--comp)
	COMP=$2
	shift
	((VAR++))
	;;
    -k|--kegg)
    KEGG=$2
    shift
    ((VAR++))
    ;;
	-t|--type)
	TYPE=$2
	shift
	((VAR++))
	;;
esac
shift
done

if [ $INPUT ] && [ $KEGG ] && [ $COMP ] && [ $TYPE ] && [ $VAR == 4 ]; then
### Check all options
if [ ! -e $INPUT ]; then 
	echo "### ERROR: no input.fofn found!"
	echo "% cp /export/EC1680U/Pipeline/Trinity/input.fofn ."
	echo "% vim input.fofn" && exit 1
elif [ ! -e $COMP ]; then
	echo "### ERROR: no comp.txt found!"
	echo "% cp /export/EC1680U/Pipeline/Trinity/comp.txt ."
	echo "% vim comp.txt" && exit 1
elif [ ! $KEGG ]; then
	echo "### ERROR: Need input KEGG organism code for KEGG annotation!"
	echo "If cannot find, please give the most nearest one"
	echo "http://www.genome.jp/kegg/catalog/org_list.html" && exit 1
elif [ ! $TYPE ]; then
    echo "### ERROR: The organism kingdom need to be euk|bac|arc !"
	echo "e.g.) -t euk | -t arc | -t bac" && exit 1
else
	echo "$(date +%T) % $0 -i $INPUT -c $COMP -k $KEGG -t $TYPE" >> "CMD_history.txt"
fi

### Check data existence
#READS=(reads/*)
#if [ ${#READS[@]} -eq 1 ]; then
#	echo "### ERROR: no read fastq existed in reads/"
#	echo "% cp /all/of/the/fastq.gz reads/" && exit 1
#fi

if [ ! -e log ]; then
	mkdir log
fi

### Run Pipeline
if [ ! -d "trinity_out_dir/Trinity.95.fasta" ]; then
	/export/EC1680U/perl/bin/Trinity/1_Trinity_assembly.pl -i $INPUT 2>&1 | tee log/1_Trinity_assembly.log
fi

if [ ! -d "1_AssemblyStats" ]; then
	/export/EC1680U/perl/bin/Trinity/2_Trinity_interval_count.pl -i $INPUT 2>&1 | tee log/2_Trinity_interval_count.log
fi

if [ ! -d "2_Quantitation" ]; then
	/export/EC1680U/perl/bin/Trinity/3_Transcript_exp_quant.pl -i $INPUT 2>&1 | tee log/3_Transcript_exp_quant.log
fi

if [ ! -d "3_DiffExpression" ]; then
	/export/EC1680U/perl/bin/Trinity/4_Diff_expression.pl -i $INPUT -c $COMP 2>&1 | tee log/4_Diff_expression.log
fi

if [ ! -d "4_Annotation" ]; then
	/export/EC1680U/perl/bin/Trinity/5_Annotation.pl -i $INPUT -c $COMP -k $KEGG -t $TYPE 2>&1 | tee log/5_Annotation.log
fi

# Report step without checking
/export/EC1680U/perl/bin/Trinity/6_Report.pl -k $KEGG 2>&1 | tee log/6_Report.log

### Help
elif [ $HELP ] || [ $VAR == 0 ]; then    
	echo "
USAGE:
### De-novo transcriptome assembly by trinity package ###
This pipeline includes 6 perl programs:
1_Trinity_assembly.pl
2_Trinity_interval_count.pl
3_Transcript_exp_quant.pl
4_Diff_expression.pl
5_Annotation.pl
6_Report.pl

[Preceding Operation]
% mkdir Trinity_NS17018-PR02
% cd Trinity_NS17018-PR02
% mkdir reads
% cp /all/of/the/fastq.gz reads/
% cp /export/EC1680U/Pipeline/Trinityinput.fofn .
% cp /export/EC1680U/Pipeline/Trinity/comp.txt .
% vim input.fofn
% vim comp.txt

% Trinity_pipe.sh -i input.fofn -c comp.txt -k sce -t euk
-----
-i|--input:
Sample information & Clean read repository path
e.g (gzip also could be used)
    ctrl_1	ctrl_1	reads/ctrl_1.R1.fastq	reads/ctrl_1.R2.fastq
	treat_1	treat_1	reads/treat_1.R1.fastq	reads/treat_1.R2.fastq
    treat_2	treat_2	reads/treat_2.R1.fastq	reads/treat_2.R2.fastq

-c|--comp:
DGE comparison file
e.g (separate with tab-delimate)
ctrl_1	treat_1
ctrl_1	treat_2
ctrl_1	treat_1,treat_2

-k|--kegg
KEGG_organism_code (http://www.genome.jp/kegg/catalog/org_list.html)

-t|--type
The organism kingdom need to be euk|bac|arc !
"
else
	echo "### ERROR: 
Parameters incorrected!" && exit 1
fi
