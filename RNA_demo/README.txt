1.) Create folder, named "reference"

2.) Download genomics / transcripts / gff3 data
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_rna.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz

3.) activate your NGS conda environment

4.) Install Trimmomatic

5.) Run the script
trimmomatic PE \
  <wt-1_R1.fastq> \
  <wt-1_R2.fastq> 
  <wt-1.R1.clean.fastq> 
  <wt-1.R1.unpaired.clean.fastq> 
  <wt-2.R2.clean.fastq> 
  <wt-2.R2.unpaired.clean.fastq> 
  ILLUMINACLIP:/export/EC1680U/software/anaconda2/pkgs/trimmomatic-0.36-5/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:1:true MINLEN:20




