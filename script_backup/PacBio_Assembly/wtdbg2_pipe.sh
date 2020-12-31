#! /bin/bash
### Goals: (subreads.fastq) ---> (contigs.polished.fasta)

#now="[`date +"%F %T"`]"

PWD=`pwd -L`

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
    fastq="$PWD/$2"
    shift
    ((VAR++))
    ;;
    -g|--genome)
    genome=$2
    shift
    ((VAR++))
    ;;
esac
shift
done

if [ $HELP ] || [ $VAR != 2 ]; then
echo "

### USAGE:
% sh wtdbg2_pipe.sh -i <fastq> -g 5m

-i/--input: subreads fastq
-g/--genome: estimated genome size
" && exit 1

fi

# Main Pipeline

echo "### Activate wtdbg2 environment..." | ts '[%Y-%m-%d %H:%M:%S]'
. activate wtdbg2

echo "### Start wtdbg2 assembly..." | ts '[%Y-%m-%d %H:%M:%S]'
if [ ! -d "wtdbg2" ]; then
	mkdir "wtdbg2"
fi

cd "wtdbg2" # change dir in

if [ ! -f "dbg.ctg.lay.gz" ]; then
	wtdbg2 -x sq -g ${genome} -i ${fastq} -t 32 -fo dbg 2>&1 | tee log.txt
fi


echo "### Generate Consensus from Layout Contigs..." | ts '[%Y-%m-%d %H:%M:%S]'
if [ ! -f "dbg.raw.ctg.fa" ]; then
	wtpoa-cns -t 32 -i dbg.ctg.lay.gz -fo dbg.raw.ctg.fa
fi

echo "### Polish Consensus..." | ts '[%Y-%m-%d %H:%M:%S]'
if [ ! -f "dbg.bam" ]; then
	minimap2 -t 32 -ax map-pb -r2k dbg.raw.ctg.fa ${fastq} | samtools sort -@32 > dbg.bam
fi
if [ ! -f "dbg.cns.ctg.fa" ]; then
	samtools view -F0x900 dbg.bam | wtpoa-cns -t 32 -d dbg.raw.ctg.fa -i - -fo dbg.cns.ctg.fa
fi

cd .. # change dir out

echo "Wtdbg2 Assembly Finish~~~ <(_ _)>" | ts '[%Y-%m-%d %H:%M:%S]'
