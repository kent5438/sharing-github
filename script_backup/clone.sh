#! /bin/bash

##### Description #####
# Use for transfering modified data to my gitlab repsitory
#######################

root_path='/Users/kentchen/Dropbox/github/NGS'
pipe_path="${root_path}/Pipeline"

### [My_script]
my_script_path="${root_path}"
server_my_script_path='/export/EC1680U/kentchen/opt/my_script';

rsync -av sheep:$server_my_script_path ${my_script_path}

### [16S]
my_16s_git_path="${root_path}/16S"
my_16s_perl_path='/export/EC1680U/perl/bin/16S'
my_16s_r_path='/export/EC1680U/Rscript/16S'
my_16s_pipe_path='/export/EC1680U/Pipeline/16S'

rsync -av --exclude='singularity' 192.168.1.221:$my_16s_perl_path/* $my_16s_git_path
rsync -av 192.168.1.221:$my_16s_r_path/* $my_16s_git_path
rsync -av 192.168.1.221:$my_16s_pipe_path $pipe_path

echo '--- 16S complete ---'

### [ITS]
my_ITS_git_path="${root_path}/ITS"
my_ITS_perl_path='/export/EC1680U/perl/bin/ITS'
my_ITS_r_path='/export/EC1680U/Rscript/ITS'

rsync -av 192.168.1.221:$my_ITS_perl_path/* $my_ITS_git_path
rsync -av 192.168.1.221:$my_ITS_r_path/* $my_ITS_git_path

### [DGE]
DGE_git_path="${root_path}/DGE"
DGE_perl_path='/export/EC1680U/perl/bin/DGE'
DGE_r_path='/export/EC1680U/Rscript/dge'
DGE_pipe_path='/export/EC1680U/Pipeline/DGE'

rsync -av 192.168.1.221:$DGE_perl_path/* $DGE_git_path
rsync -av 192.168.1.221:$DGE_r_path/* $DGE_git_path
rsync -av  192.168.1.221:$DGE_pipe_path/multiqc_config.yaml $pipe_path/DGE

echo '--- DGE complete ---'

### [Trinity]
trinity_git_path="${root_path}/Trinity"
trinity_perl_path='/export/EC1680U/perl/bin/Trinity'
trinity_r_path='/export/EC1680U/Rscript/Trinity'
trinity_pipe_path='/export/EC1680U/Pipeline/Trinity'

rsync -av 192.168.1.221:$trinity_perl_path/* $trinity_git_path
rsync -av 192.168.1.221:$trinity_r_path/* $trinity_git_path/Rscript
rsync -av 192.168.1.221:$trinity_pipe_path $pipe_path

echo '--- Trinity complete ---'

### [PacBio_Assembly_Denovo_Assembly]
pb_asm_git_path="${root_path}/PacBio_Assembly"
pb_asm_perl_path='/export/EC1680U/perl/bin/Pacbio'
pb_asm_r_path='/export/EC1680U/Rscript/Pacbio'
pb_asm_pipe_path='/export/EC1680U/Pipeline/Pacbio_Denovo/'

rsync -av 192.168.1.221:$pb_asm_perl_path/* $pb_asm_git_path
rsync -av 192.168.1.221:$pb_asm_r_path/* $pb_asm_git_path/Rscript
rsync -av 192.168.1.221:$pb_asm_pipe_path $pipe_path

echo '--- PB denovo assembly complete ---'

### [MiSeq_Denovo_Assembly]
miseq_asm_git_path="${root_path}/Illumina_Assembly"
miseq_asm_perl_path='/export/EC1680U/perl/bin/Assembly'
miseq_asm_pipe_path='/export/EC1680U/Pipeline/Miseq_Assembly'

rsync -av 192.168.1.221:$miseq_asm_perl_path/* $miseq_asm_git_path
rsync -av 192.168.1.221:$miseq_asm_pipe_path $pipe_path

echo '--- Miseq denovo assembly complete ---'

### [10x_supernova]
my_10x_supernova_git_path="${root_path}/10x_Supernova"
my_10x_supernova_perl_path='/export/EC1680U/perl/bin/10x_Supernova'
my_10x_supernova_pipe_path='/export/EC1680U/Pipeline/10x_Supernova'

rsync -av 192.168.1.221:$my_10x_supernova_perl_path/* $my_10x_supernova_git_path
rsync -av 192.168.1.221:$my_10x_supernova_pipe_path $pipe_path

echo '--- 10x supernova complete ---'

### [Amplicon]
amp_git_path="${root_path}/Amplicon/"
amp_perl_path='/export/EC1680U/perl/bin/Amplicon'
amp_r_path='/export/EC1680U/Rscript/Amplicon'

rsync -av 192.168.1.221:$amp_perl_path/*.pl $amp_git_path
#rsync -avq 192.168.1.221:$amp_perl_path/*.sh $amp_git_path
#rsync -avq 192.168.1.221:$amp_perl_path/*.R $amp_git_path/Rscript

echo ''
echo '--- All Complete ---'
echo ''
