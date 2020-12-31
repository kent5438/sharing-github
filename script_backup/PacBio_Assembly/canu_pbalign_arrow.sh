#! /bin/sh

JOB_ROOT='/export/EC2480U/smrtlink_userdata_Ox/jobs_root/001/'
VAR=0

while [ $# -gt 0 ]
do
OPT=$1

case $OPT in
    -h|--help)
    HELP=$1
    ((VAR++))
    ;;
    -i|--id)
    SMRT_ID=$2
    shift
    ((VAR++))
    ;;
    -p|--prefix)
    PREFIX=$2
    shift
    ((VAR++))
    ;;
    -g|--genome)
    GENOME=$2
    shift
    ((VAR++))
    ;;
esac
shift
done

if [ $HELP ] || [ $VAR != 3 ]; then
echo "

### USAGE:
% . activate pb-rsii-denovo
% sh canu_pbalign_arrow.sh -i 001701 -p S-58 -g 2.35m

" && exit 1

fi

### Canu pipe
echo "#!/bin/bash
#SBATCH --job-name=canu
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8
#SBATCH --output=canu.log
#SBATCH --error=canu.err

canu \
-p ${PREFIX} -d canu \
genomeSize=${GENOME} \
-pacbio-raw ${JOB_ROOT}/${SMRT_ID}/tasks/pbcoretools.tasks.bam2fastq_archive-0/reads.fastq \
useGrid=true gridEngine="slurm"
" > 'canu_slurm.sh'
if [ ! -e "canu/${PREFIX}.contigs.fasta" ]; then
	srun canu_slurm.sh
fi

cd canu/

### pbalign pipe
echo "#!/bin/bash
#SBATCH --job-name=pbalign
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8
#SBATCH --output=pbalign.log
#SBATCH --error=pbalign.err

/opt/smrtlink/smrtcmds/bin/pbalign \
--nproc 8 --tmpDir ./tmp \
-v ${JOB_ROOT}/${SMRT_ID}/tasks/pbcoretools.tasks.filterdataset-0/filtered.subreadset.xml \
${PREFIX}.contigs.fasta ${PREFIX}.contigs.align.bam
" > 'pbalign_slurm.sh'

if [ ! -e "${PREFIX}.contigs.align.bam" ]; then
	srun pbalign_slurm.sh
fi

### Arrow pipe
echo "#!/bin/bash
#SBATCH --job-name=arrow
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8
#SBATCH --output=arrow.log
#SBATCH --error=arrow.err

/opt/smrtlink/smrtcmds/bin/arrow \
${PREFIX}.contigs.align.bam \
-r ${PREFIX}.contigs.fasta \
-o ${PREFIX}.var.gff \
-o ${PREFIX}.polished.fasta \
-o ${PREFIX}.polished.fastq \
-j 8
" > 'arrow_slurm.sh'


if [ ! -e "${PREFIX}.polished.fasta" ]; then
	srun arrow_slurm.sh
fi

cd ..
