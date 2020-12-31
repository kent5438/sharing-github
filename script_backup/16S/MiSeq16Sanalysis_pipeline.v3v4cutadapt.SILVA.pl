#! /usr/bin/perl -w

# Written by Keh-Ming Wu, Ver 1.01, 2014.10.15
#
# 16S Pipeline after QC step of sequencing


use Getopt::Long;
#use Getopt::Std;
#use strict;
use FileHandle;
use Spreadsheet::WriteExcel;
use PDL;

STDOUT->autoflush(1);

my $qiime_docker = "sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` bwawrik/qiime:latest";

##===================================================
MainPrg:

$Usage =	
    "=======================================\n".
		"* Something Wrong with the parameters *\n".
		"=======================================\n\n".
		"Usage:\n".
		"$0 -s sampleNameList.txt -p [primerSeq.txt] -r [readStatTag:Y/N, default:Y] -o [OTUstatTag:Y/N, default:N] -i [baseSpaceTag:Y/N,default:N]\n".
		"$0 -s analysisTarget20.txt -p ../V1V3primers.alienSeq -n N\n".
		"$0 -s MS14020_ST_44samples_list.txt -p V1V3primers.alienSeq -r N -o Y -i Y\n".
		"$0 -s MS14021_SI_38samples_list.txt -p V1V3primers.alienSeq -r N -o Y\n".
		"$0 -s MS14005_6samples_target1.txt -p 515F-884R_primers.alienSeq\n".
		"$0 -s MS14005_6samples_target2.txt -p 789F-1053R_primers.alienSeq\n".
		"$0 -s MS14005_6samples_target3.txt -p 926F-1392R_primers.alienSeq\n".
		"$0 -s MS14022_8samples_target1.txt -p 515F-884R_primers.alienSeq\n".
		"$0 -s MS14022_8samples_target2.txt -p 789F-1053R_primers.alienSeq\n".
		"$0 -s MS14022_8samples_target3.txt -p 926F-1392R_primers.alienSeq\n\n";

#default subfolder for each project
$folder[0]="0_cleanFastq";
$folder[1]="1_joinPairs";
$folder[2]="2_removePrimers";
$folder[3]="3_screenQuality";
$folder[4]="4_checkChimera";
$folder[5]="5_assignTaxa";
$folder[6]="6_biom";
$folder[7]="meta16S_OUTPUT";
$OTUdataFolder="$folder[7]/OTUs_eachSample";
$folder[8]="meta16S_baseSpaceSourceFastq";


if (!-d $folder[0])
{
	print "Please put clean fastq files into '$folder[0]', then execute this script again!~\n\n"; 
	exit;
}



GetOptions(
        's=s'  => \$sFile,
        'p=s'  => \$primerSet,
        'r=s'  => \$readStatTag,
        'o=s'  => \$OTUstatTag,
        'i=s'  => \$baseSpaceTag
);

open OUT, ">>CMD_history.txt";
print OUT getCurrentDateTime(). "$0 -s $sFile -p $primerSet -r $readStatTag -o $OTUstatTag -i $baseSpaceTag 2>&1 | tee mothur_pipe.log\n";
close OUT;


die "$Usage" if (!defined($sFile) && !defined($primerSet));

if (defined($sFile))
{
	if (!-e $sFile){die "$sFile not found!!\n\n";}
}
else
{
	die "sampleNameList.txt not provided with -s parameter! \n\n";	
}

if (defined($primerSet))
{
	if (!-e $primerSet){die "$primerSet not found!!\n\n";}
}
else
{
	if (-e "../V1V3primers.alienSeq")
	{
		print "primerSeq.txt NOT provided, but found '../V1V3primers.alienSeq'! use it instead!!\n\n";
	}
	else
	{
		die "primerSeq.txt not provided with -p parameter! \n\n";	
	}
}

if (defined($readStatTag))
{
	if ($readStatTag !~ /^[YN]/i ){die "readStatTag should be either Y or N only!!\n\n";}
	elsif ($readStatTag =~ /^[Y]/i){ $readStatTag="Y"; }
	else { $readStatTag="N"; }
}
else
{
	$readStatTag="Y";
	print "Should perform read info statistics later...\n\n";
}

if (defined($OTUstatTag))
{
	if ($OTUstatTag !~ /^[YN]/i ){die "OTUstatTag should be either Y or N only!!\n\n";}
	elsif ($OTUstatTag =~ /^[Y]/i){ $OTUstatTag="Y"; }
	else { $OTUstatTag="N"; }
}
else
{
	$OTUstatTag="N";
	print "Should perform read info statistics later...\n\n";
}

if (defined($baseSpaceTag))
{
	if ($baseSpaceTag !~ /^[YN]/i ){die "baseSpaceTag should be either Y or N only!!\n\n";}
	elsif ($baseSpaceTag =~ /^[Y]/i){ $baseSpaceTag="Y"; }
	else { $baseSpaceTag="N"; }
}
else
{
	$baseSpaceTag="N";
	#print "Should prepare fastq.gz for BaseSpace 16S APP. later...\n\n";
}


foreach $folder (@folder)
{
	if (!-d $folder){print "Folder $folder for the 16S analysis is not found! Make folder now~\n"; mkdir ($folder);}
	
	if ($folder eq $folder[7]){if (!-d $OTUdataFolder){mkdir ($OTUdataFolder);}}
}




$ref2ChimeraUchime="/export/EC1680U/metagenomics/16Sdb/SILVA/silva.bacteria/silva.gold.ng.fasta"; # 5181
if (!-e $ref2ChimeraUchime){die "$ref2ChimeraUchime not found!!\n\n"}


#$refFastaDB[0]="/export/EC1680U/metagenomics/16Sdb/RDP/trainset9_032012.pds.fasta";
#$refTaxonomy[0]="/export/EC1680U/metagenomics/16Sdb/RDP/trainset9_032012.pds.tax";

# "/export/EC1680U/metagenomics/16Sdb/MothurUsage/greenGenes/gg_13_8_99.refalign"
#$refAlignDB[0]="/export/EC1680U/metagenomics/16Sdb/SILVA/silva_180117/silva.raw_taxonomy.align"; # 213,119
$refAlignDB[1]="/export/EC1680U/DataBase/16S/SILVA_132_release/rep_set_aligned/99/99_alignment.fna"; # 203,452 
if (!-e $refAlignDB[1]){die "$refAlignDB[1] not found!!\n\n"}

$refFastaDB[1]="/export/EC1680U/DataBase/16S/SILVA_132_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna";
#$refFastaDB[2]="/export/EC1680U/metagenomics/16Sdb/MothurUsage/gg_13_8_99.fasta"; #
if (!-e $refFastaDB[1]){die "$refFastaDB[1] not found!!\n\n"}

$refTaxonomy[1]="/export/EC1680U/DataBase/16S/SILVA_132_release/taxonomy/16S_only/99/raw_taxonomy.txt"; # 4	Root;k__Archaea;p__[Parvarchaeota];c__[Parvarchaea];o__YLA114;
#$refTaxonomy[2]="/export/EC1680U/metagenomics/16Sdb/MothurUsage/gg_13_8_99.gg.tax"; # 4	k__Archaea;p__[Parvarchaeota];c__[Parvarchaea];o__YLA114;
if (!-e $refTaxonomy[1]){die "$refTaxonomy[1] not found!!\n\n"}

my $silva_ref_tree="/export/EC1680U/DataBase/16S/SILVA_132_release/trees/99/99_otus.tre";

#$refotuMapTable[1]="/export/EC1680U/metagenomics/16Sdb/MothurUsage/greenGenes/GG_13_5_otuMapTable/99_otu_map.txt";
#if (!-e $refotuMapTable[1]){die "$refotuMapTable[1] not found!!\n\n"}

#====================
$sampleIDerr=0;$sCount=0;

print getCurrentDateTime()."\t";
print "Reading sampleList File.... $sFile\n";
open(IN, $sFile) || die "The file $sFile can not be openned.\n";
while(defined($line=<IN>))
{ 
	chomp $line;  $line=~ s/\r\n?/\n/g;
	next if ($line =~ /^\#/);
	next if ($line eq "");
	
	@sArray=split /\t/,$line;
	#$itemNum=$#sArray +1; 
	$sampleID=trim($sArray[0]);
	#$bcArray[0]=trim($sArray[0]);
	
	$sCount++;
	print "$sCount|$sampleID|\n";
	if ($sampleID !~ /[A-Za-z0-9\-_]/){print "!! SampleID '$sampleID' should contain only alphanumeric characters & dash(-), or underline(_)! "; $sampleIDerr=1;}
	
	
	$sampleWithDash=$sampleID;
	if ($sampleWithDash =~ /-/){$sampleWithDash=~ s/-/\\-/g; }
	
	push @sampleIDs,$sampleID;
	push @sampleIDDash,$sampleWithDash;
	
	
}
close (IN);

if ($sampleIDerr==1) {die "\n\nAt least one sampleID is invalid!! Please check the sampleList!~\n\n";}



$fastqSuffix[0]="clean.fastq";
$fastqSuffix[1]="assembled.fastq";
$fastqSuffix[2]="assembledTrimPrimers.fastq";

$fastqSuffix[3]="qualified.fastq";
$fastaSuffix[3]="qualified.fasta";
$qualSuffix[3]="qualified.qual";

$fastaSuffix[4]="effective.fasta";
$qualSuffix[4]="effective.qual";

$folderSuffix[5]="16Sanalysis";
$fastaSuffix[5]="fasta";
$qualSuffix[5]="qual";


#$fastaSuffix="effective.fasta";
#$qualSuffix="effective.qual";


# check files if exist
$fileExistCK=0;
print "\n";
print getCurrentDateTime()."\t";
print "Checking file existence..\n";
for ($i=0;$i< $sCount ;$i++)
{
	# Step0: prepare the source fastq (clean.fastq)
	if (-e "./$folder[0]/$sampleIDs[$i].R1.$fastqSuffix[0]" && -e "./$folder[0]/$sampleIDs[$i].R2.$fastqSuffix[0]")
	{
		$cleanPEgzTag=0;
		$fileExistCK=0;
	}
	elsif (-e "./$folder[0]/$sampleIDs[$i].R1.$fastqSuffix[0].gz" && -e "./$folder[0]/$sampleIDs[$i].R2.$fastqSuffix[0].gz")
	{
		$cleanPEgzTag=1;
		$fileExistCK=0;
	}
	else 
	{
		print "./$folder[0]/$sampleIDs[$i].R1.$fastqSuffix[0](or .gz) or ./$folder[0]/$sampleIDs[$i].R2.$fastqSuffix[0](or .gz) not found!!\n";
		$fileExistCK=1;
		exit 1;
	}

	# Step1: prepare joined reads (assembled.fastq)
	if (-e "./$folder[1]/$sampleIDs[$i].$fastqSuffix[1]" )
	{
		$fileExistCK=0;
		$joinedGzTag=0;
	}
	elsif (-e "./$folder[1]/$sampleIDs[$i].$fastqSuffix[1].gz" )
	{
		$fileExistCK=0;
		$joinedGzTag=1;
	}
	else 
	{

		if ($cleanPEgzTag==1)
		{
			`gzip -dc ./$folder[0]/$sampleIDs[$i].R1.$fastqSuffix[0].gz > ./$folder[0]/$sampleIDs[$i].R1.$fastqSuffix[0]`;
			`gzip -dc ./$folder[0]/$sampleIDs[$i].R2.$fastqSuffix[0].gz > ./$folder[0]/$sampleIDs[$i].R2.$fastqSuffix[0]`;
		}
		chdir ("./$folder[1]");
		
		#old: pear -m 600 -n 300 -u 0 -v 15 -j 24 -o $sampleIDs[$i] -f
		#minimal 100 bp for joined fragments
		print "pear -m 500 -n 50 -s 2 -u 1 -p 0.05 -v 15 -j 24 -o $sampleIDs[$i] -f ../$folder[0]/$sampleIDs[$i].R1.$fastqSuffix[0] -r ../$folder[0]/$sampleIDs[$i].R2.$fastqSuffix[0] &> $sampleIDs[$i].pear.log\n";
	  	    `/export/EC1680U/software/pear-0.9.10-bin-64/pear -m 500 -n 50 -s 2 -u 1 -p 0.05 -v 15 -j 24 -o $sampleIDs[$i] -f ../$folder[0]/$sampleIDs[$i].R1.$fastqSuffix[0] -r ../$folder[0]/$sampleIDs[$i].R2.$fastqSuffix[0] &> $sampleIDs[$i].pear.log`;

		chdir ("..");
		if ($cleanPEgzTag==1)
		{
			# if gz file exists, then delete the uncompressed fastq file
			if (-e "./$folder[0]/$sampleIDs[$i].R1.$fastqSuffix[0].gz")
			{
				unlink ("./$folder[0]/$sampleIDs[$i].R1.$fastqSuffix[0]");
			}
			if (-e "./$folder[0]/$sampleIDs[$i].R2.$fastqSuffix[0].gz")
			{
				unlink ("./$folder[0]/$sampleIDs[$i].R2.$fastqSuffix[0]");
			}
		}
		
		if (-e "./$folder[1]/$sampleIDs[$i].$fastqSuffix[1]" )
		{
			$fileExistCK=0;
			$joinedGzTag=0;
		}
		elsif (-e "./$folder[1]/$sampleIDs[$i].$fastqSuffix[1].gz" )
		{
			$fileExistCK=0;
			$joinedGzTag=1;
		}
		else
		{
			print "./$folder[1]/$sampleIDs[$i].$fastqSuffix[1](or .gz) not found!!\n";
			$fileExistCK=1;
			exit 1;
		}
		
	
	}
	
	# Step2: prepare primer-trimmed reads (assembledTrimPrimers.fastq)
	if (-e "./$folder[2]/$sampleIDs[$i].$fastqSuffix[2]" )
	{
		$fileExistCK=0;
		$assembledTrimPrimersGzTag=0;
	}
	elsif (-e "./$folder[2]/$sampleIDs[$i].$fastqSuffix[2].gz" )
	{
		$fileExistCK=0;
		$assembledTrimPrimersGzTag=1;
	}
	else 
	{
		if ($joinedGzTag==1)
		{
			`gzip -dc ./$folder[1]/$sampleIDs[$i].$fastqSuffix[1].gz > ./$folder[1]/$sampleIDs[$i].$fastqSuffix[1]`;
		}
		chdir ("./$folder[2]");
		
		
		$trimOption=1;
		if ($trimOption== 1)
		{
			#V1V3primers.alienSeq 
			#AGAGTTTGATCMTGGCTCAG
			#GTATTACCGCGGCKGCTG  => RC => CAGCMGCCGCGGTAATAC
			
			if ($primerSet =~ /V3V4/)
			{
				`/export/EC1680U/software/anaconda2/bin/cutadapt -g CCTACGGGNGGCWGCAG -a GACTACHVGGGTATCTAATCC -m 200 --discard-untrimmed -o $sampleIDs[$i].$fastqSuffix[2] ../$folder[1]/$sampleIDs[$i].$fastqSuffix[1] &> $sampleIDs[$i].assembledTrimPrimers.log`;
#				print "AlienTrimmer -q 0 -l 300 -p 50 -i ../$folder[1]/$sampleIDs[$i].$fastqSuffix[1] -o $sampleIDs[$i].$fastqSuffix[2] -c ../$primerSet &> $sampleIDs[$i].assembledTrimPrimers.log\n";
#				`/export/EC1680U/software/AlienTrimmer_0.4.0/src/AlienTrimmer -q 0 -l 300 -p 50 -i ../$folder[1]/$sampleIDs[$i].$fastqSuffix[1] -o $sampleIDs[$i].$fastqSuffix[2] -c ../$primerSet &> $sampleIDs[$i].assembledTrimPrimers.log`;
				`cat $sampleIDs[$i].assembledTrimPrimers.log`;
			}
			else
			{
				print "AlienTrimmer -q 0 -l 200 -p 80 -i ../$folder[1]/$sampleIDs[$i].$fastqSuffix[1] -o $sampleIDs[$i].$fastqSuffix[2] -c ../$primerSet &> $sampleIDs[$i].assembledTrimPrimers.log\n";
				`/export/EC1680U/software/AlienTrimmer_0.4.0/src/AlienTrimmer -q 0 -l 200 -p 80 -i ../$folder[1]/$sampleIDs[$i].$fastqSuffix[1] -o $sampleIDs[$i].$fastqSuffix[2] -c ../$primerSet &> $sampleIDs[$i].assembledTrimPrimers.log`;
				`cat $sampleIDs[$i].assembledTrimPrimers.log`;
			}
		}
		elsif ($trimOption== 2)
		{
		#cutadapt -a ADAPTER_FWD -o trimmed.1.fastq reads1.fastq
		#cutadapt -a ADAPTER_REV -o trimmed.2.fastq reads2.fastq
		#cutadapt -g AGAGTTTGATCMTGGCTCAG --match-read-wildcards --discard-untrimmed -O 15 -o sample014.trimP1.fastq  sample014.assembled.fastq
		#
		#O must reverse&complement the 3'primerSeq.  
		#cutadapt -a CAGCMGCCGCGGTAATAC --match-read-wildcards --discard-untrimmed -O 15 -m 200 -o sample014.trimmed.fastq  sample014.trimP1.fastq 
		
		#Not good at trimming 3' end
		#cutadapt -g AGAGTTTGATCMTGGCTCAG -a CAGCMGCCGCGGTAATAC  --match-read-wildcards --discard-untrimmed -O 15 -o sample014.trimP1.fastq  sample014.assembled.fastq
		
		}
		
		chdir ("..");
		if ($joinedGzTag==1)
		{
			unlink ("./$folder[1]/$sampleIDs[$i].$fastqSuffix[1]");
		}
		
		if (-e "./$folder[2]/$sampleIDs[$i].$fastqSuffix[2]" )
		{
			$fileExistCK=0;
			$assembledTrimPrimersGzTag=0;
		}
		elsif (-e "./$folder[2]/$sampleIDs[$i].$fastqSuffix[2].gz" )
		{
			$fileExistCK=0;
			$assembledTrimPrimersGzTag=1;
		}
		else
		{
			print "./$folder[2]/$sampleIDs[$i].$fastqSuffix[2](or .gz) not found!!\n";
			$fileExistCK=1;
			exit 1;
		}

	}
	
	#step3: prepare qualified fastq (qualified.fastq) and fasta+qual (qualified.fasta & qualified.qual)
	if (-e "./$folder[3]/$sampleIDs[$i].$fastqSuffix[3]" )
	{
		$fileExistCK=0;
		$qualifiedFQGzTag=0;
	}
	elsif (-e "./$folder[3]/$sampleIDs[$i].$fastqSuffix[3].gz" )
	{
		$fileExistCK=0;
		$qualifiedFQGzTag=1;
	}
	else 
	{
		
		if ($assembledTrimPrimersGzTag==1)
		{
			`gzip -dc ./$folder[2]/$sampleIDs[$i].$fastqSuffix[2].gz > ./$folder[2]/$sampleIDs[$i].$fastqSuffix[2]`;
		}
		chdir ("./$folder[3]");
		
		#Step3-1: Quality trimming & length filtering
		if ($primerSet =~ /V3V4/)
		{
			print "java -jar /export/EC1680U/software/trimmomatic-0.33.jar SE -threads 24 -phred33 -trimlog temp_$sampleIDs[$i].log ../$folder[2]/$sampleIDs[$i].$fastqSuffix[2] $sampleIDs[$i].$fastqSuffix[3] SLIDINGWINDOW:10:10 MINLEN:300\n";		
		  	    `java -jar /export/EC1680U/software/trimmomatic-0.33.jar SE -threads 24 -phred33 -trimlog temp_$sampleIDs[$i].log ../$folder[2]/$sampleIDs[$i].$fastqSuffix[2] $sampleIDs[$i].$fastqSuffix[3] SLIDINGWINDOW:10:10 MINLEN:300`;
		}
		else
		{
			print "java -jar /export/tools/trimmomatic.jar SE -threads 24 -phred33 -trimlog temp_$sampleIDs[$i].log ../$folder[2]/$sampleIDs[$i].$fastqSuffix[2] $sampleIDs[$i].$fastqSuffix[3] SLIDINGWINDOW:10:10 MINLEN:200\n";		
		  	    `java -jar /export/tools/trimmomatic.jar SE -threads 24 -phred33 -trimlog temp_$sampleIDs[$i].log ../$folder[2]/$sampleIDs[$i].$fastqSuffix[2] $sampleIDs[$i].$fastqSuffix[3] SLIDINGWINDOW:10:10 MINLEN:200`;
		}
#		system('perl /export/EC1680U/jsyf_16S/fastqFilter.pl '.$sampleIDs[$i].'.'.$fastqSuffix[3].' 27');
		#Step3-2: Convertion of fasta/qual from fastq
		print "mothur \"#fastq.info(format=illumina1.8+,qfile=T,fastq=$sampleIDs[$i].$fastqSuffix[3])\"\n";
					`/export/EC1680U/software/anaconda2/bin/mothur "#fastq.info(format=illumina1.8+,qfile=T,fastq=$sampleIDs[$i].$fastqSuffix[3])"`;

		chdir ("..");
		if ($assembledTrimPrimersGzTag==1)
		{
			unlink ("./$folder[2]/$sampleIDs[$i].$fastqSuffix[2]");
			unlink ("./$folder[3]/temp_$sampleIDs[$i].log");
		}
		
		if (-e "./$folder[3]/$sampleIDs[$i].$fastqSuffix[3]" )
		{
			$fileExistCK=0;
			$qualifiedFQGzTag=0;
		}
		elsif (-e "./$folder[3]/$sampleIDs[$i].$fastqSuffix[3].gz" )
		{
			$fileExistCK=0;
			$qualifiedFQGzTag=1;
		}
		else
		{
			print "./$folder[3]/$sampleIDs[$i].$fastqSuffix[3](or .gz) not found!!\n";
			$fileExistCK=1;
			exit 1;
		}
	}
	
	#Step3-3: check qualified FASTA
	if (-e "./$folder[3]/$sampleIDs[$i].$fastaSuffix[3]" )
	{
		$fileExistCK=0;
		$qualifiedFAGzTag=0;
	}
	elsif (-e "./$folder[3]/$sampleIDs[$i].$fastaSuffix[3].gz" )
	{
		$fileExistCK=0;
		$qualifiedFAGzTag=1;
	}
	else 
	{
		print "./$folder[3]/$sampleIDs[$i].$fastaSuffix[3](or .gz) not found!!!\n";
		$fileExistCK=1;
		
		if (-e "./$folder[3]/$sampleIDs[$i].$fastqSuffix[3]")
		{
			print "Try to delete ./$folder[3]/$sampleIDs[$i].$fastqSuffix[3] and run this script again!!\n\n";
		}
		exit 1;
	}
	#Step3-3: check qualified QUAL
	if (-e "./$folder[3]/$sampleIDs[$i].$qualSuffix[3]" )
	{
		$fileExistCK=0;
		$qualifiedQAGzTag=0;
	}
	elsif (-e "./$folder[3]/$sampleIDs[$i].$qualSuffix[3].gz" )
	{
		$fileExistCK=0;
		$qualifiedQAGzTag=1;
	}
	else 
	{
		print "./$folder[3]/$sampleIDs[$i].$qualSuffix[3](or .gz) not found!!\n";
		$fileExistCK=1;
		
		if (-e "./$folder[3]/$sampleIDs[$i].$fastqSuffix[3]")
		{
			print "Try to delete ./$folder[3]/$sampleIDs[$i].$fastqSuffix[3] and run this script again!!\n\n";
		}
		exit 1;
	}



	#Step4: Chimera check and filtering
	if (!-e "./$folder[3]/$sampleIDs[$i].qualified.uchime.accnos" || !-e "./$folder[3]/$sampleIDs[$i].qualified.uchime.chimeras")
	{
		print "./$folder[3]/$sampleIDs[$i].qualified.uchime.accnos or ./$folder[3]/$sampleIDs[$i].qualified.uchime.chimeras not found!!\n";
		#$fileExistCK=1;
	
		print "\nPerform Mothur 'chimera.uchime' step automatically... \n\n";
			
			chdir ("./$folder[3]");
			print "mothur \"#chimera.uchime(processors=30,fasta=$sampleIDDash[$i].$fastaSuffix[3], reference=$ref2ChimeraUchime)\"\n";
			print "mothur \"#remove.seqs(accnos=$sampleIDs[$i].qualified.uchime.accnos,fasta=$sampleIDs[$i].$fastaSuffix[3],qfile=$sampleIDs[$i].$qualSuffix[3])\"\n";
			`/export/EC1680U/software/anaconda2/bin/mothur "#chimera.uchime(processors=30,fasta=$sampleIDDash[$i].$fastaSuffix[3], reference=$ref2ChimeraUchime)" &> /dev/null`;
			`mv $sampleIDs[$i].qualified.ref.uchime.accnos $sampleIDs[$i].qualified.uchime.accnos`;
			`mv $sampleIDs[$i].qualified.ref.uchime.chimeras $sampleIDs[$i].qualified.uchime.chimeras`;
			`/export/EC1680U/software/anaconda2/bin/mothur "#remove.seqs(accnos=$sampleIDs[$i].qualified.uchime.accnos,fasta=$sampleIDs[$i].$fastaSuffix[3],qfile=$sampleIDs[$i].$qualSuffix[3])"`;
			`mv $sampleIDs[$i].qualified.pick.fasta ../$folder[4]/$sampleIDs[$i].$fastaSuffix[4]`;
			`mv $sampleIDs[$i].qualified.pick.qual ../$folder[4]/$sampleIDs[$i].$qualSuffix[4]`;
			`rm -rf $sampleIDs[$i].qualified.fasta*.uchime_formatted`;
			`rm -rf $sampleIDs[$i].qualified.uchime*.temp`;
			chdir ("..");

	}

	if (-e "./$folder[4]/$sampleIDs[$i].$fastaSuffix[4]" )
	{
		$fileExistCK=0;
		$effectiveFAGzTag=0;
	}
	elsif (-e "./$folder[4]/$sampleIDs[$i].$fastaSuffix[4].gz" )
	{
		$fileExistCK=0;
		$effectiveFAGzTag=1;
	}
	else 
	{
		print "./$folder[4]/$sampleIDs[$i].$fastaSuffix[4](or .gz) not found!!\n";
		$fileExistCK=1;
		exit 1;
	}
	if (-e "./$folder[4]/$sampleIDs[$i].$qualSuffix[4]" )
	{
		$fileExistCK=0;
		$effectiveQAGzTag=0;
	}
	elsif (-e "./$folder[4]/$sampleIDs[$i].$qualSuffix[4].gz" )
	{
		$fileExistCK=0;
		$effectiveQAGzTag=1;
	}
	else 
	{
		print "./$folder[4]/$sampleIDs[$i].$qualSuffix[4](or .gz) not found!!\n";
		$fileExistCK=1;
		exit 1;
	}
	
	
	#if (!-e "./$folder[4]/$sampleIDs[$i].$fastaSuffix[4]" || !-e "./$folder[4]/$sampleIDs[$i].$qualSuffix[4]")
	#{
	#	print "./$folder[4]/$sampleIDs[$i].$fastaSuffix[4] or ./$folder[4]/$sampleIDs[$i].$qualSuffix[4] not found!!\n";
	#	
	#	$fileExistCK=1;
	#	exit 1;
	#}
	
}
print "\n";
#if ($fileExistCK==1){print "\n**At least one file not found**\nPress Y to CONTINUE, others to STOP...\n\n";$ifCS=<STDIN>; if ($ifCS !~ /Y/i){exit;}}


if ($readStatTag eq "Y")
{
	print "\nStat. the number and quality of files at each folder...\n\n";
	print getCurrentDateTime()."\t";
	print "Analysing stat infos..\n";
	# count reads and stats
	for ($f=0;$f< 5 ;$f++)
	{
		print "\nProcessing $folder[$f] ...\n";
		
		for ($i=0;$i< $sCount ;$i++)
		{
			print ($i+1);
			print "\t$sampleIDs[$i] ~\n";
			
			if ($f==0) #clean.fastq
			{
				if ($cleanPEgzTag==1)
				{($total,$totalLen,$avgLen)=statFastq("$folder[0]/$sampleIDs[$i].R1.$fastqSuffix[0].gz",$cleanPEgzTag);}
				else
				{($total,$totalLen,$avgLen)=statFastq("$folder[0]/$sampleIDs[$i].R1.$fastqSuffix[0]",$cleanPEgzTag);}
				$cleanReadCount{$sampleIDs[$i]}=$total;
				$cleanBaseCount{$sampleIDs[$i]}=$totalLen;
				
				print "$sampleIDs[$i].R1.$fastqSuffix[0]: total $cleanReadCount{$sampleIDs[$i]} reads, $cleanBaseCount{$sampleIDs[$i]} bases, average $avgLen bp in length\n";
			}
			if ($f==1) #joined fastq
			{
				if ($joinedGzTag==1)
				{($total,$totalLen,$avgLen)=statFastq("$folder[1]/$sampleIDs[$i].$fastqSuffix[1].gz",1);}
				else
				{($total,$totalLen,$avgLen)=statFastq("$folder[1]/$sampleIDs[$i].$fastqSuffix[1]",0);}
				
				$joinedReadCount{$sampleIDs[$i]}=$total;
				$joinedBaseCount{$sampleIDs[$i]}=$totalLen;
				
				print "$sampleIDs[$i].$fastqSuffix[1]: total $joinedReadCount{$sampleIDs[$i]} reads, $joinedBaseCount{$sampleIDs[$i]} bases, average $avgLen bp in length\n";
			}		
			if ($f==3) # qualified fastq
			{
				if ($qualifiedFQGzTag==1)
				{($total,$totalLen,$avgLen)=statFastq("$folder[3]/$sampleIDs[$i].$fastqSuffix[3].gz",1);}
				else
				{($total,$totalLen,$avgLen)=statFastq("$folder[3]/$sampleIDs[$i].$fastqSuffix[3]",0);}
				
				$qualifiedReadCount{$sampleIDs[$i]}=$total;
				$qualifiedBaseCount{$sampleIDs[$i]}=$totalLen;
				
				print "$sampleIDs[$i].$fastqSuffix[3]: total $qualifiedReadCount{$sampleIDs[$i]} reads, $qualifiedBaseCount{$sampleIDs[$i]} bases, average $avgLen bp in length\n";
			}
			
			if ($f==4) # non-chimera fasta & Q20 Q30 GC stats
			{
				if ($effectiveFAGzTag==1)
				{($total,$totalLen,$avgLen,$GCpct,$mean,$mode,$median,$min,$max,$stdev)=statFasta("$folder[4]/$sampleIDs[$i].$fastaSuffix[4].gz",1);}
				else
				{($total,$totalLen,$avgLen,$GCpct,$mean,$mode,$median,$min,$max,$stdev)=statFasta("$folder[4]/$sampleIDs[$i].$fastaSuffix[4]");}
				
				if ($effectiveQAGzTag==1)
				{($avgQ, $Q20pct,$Q30pct)=statQual("$folder[4]/$sampleIDs[$i].$qualSuffix[4].gz",1);}
				else
				{($avgQ, $Q20pct,$Q30pct)=statQual("$folder[4]/$sampleIDs[$i].$qualSuffix[4]");}
				
				$effectiveReadCount{$sampleIDs[$i]}=$total;
				$effectiveBaseCount{$sampleIDs[$i]}=$totalLen;
				$effectiveAvgLen{$sampleIDs[$i]}=$avgLen;
				$effectiveGCpct{$sampleIDs[$i]}=$GCpct;
				$effectiveAvgQ{$sampleIDs[$i]}=$avgQ;
				$effectiveQ20pct{$sampleIDs[$i]}=$Q20pct;
				$effectiveQ30pct{$sampleIDs[$i]}=$Q30pct;
				
				$effectiveMeanLen{$sampleIDs[$i]}=$mean;
				$effectiveModeLen{$sampleIDs[$i]}=$mode;
				$effectiveMedianLen{$sampleIDs[$i]}=$median;
				$effectiveMinLen{$sampleIDs[$i]}=$min;
				$effectiveMaxLen{$sampleIDs[$i]}=$max;
				$effectiveStdevLen{$sampleIDs[$i]}=$stdev;

				$effectivePct{$sampleIDs[$i]}=sprintf("%.02f",$total/$cleanReadCount{$sampleIDs[$i]}*100);
				
				print "$sampleIDs[$i].$fastaSuffix[4]: total $total reads, $totalLen bases, average $avgLen bp in length\n";
				print "$sampleIDs[$i].$fastaSuffix[4]: average Q value= $avgQ, Q20%= $Q20pct, Q30%= $Q30pct, GC%=$GCpct\n";
			}
		}
	}
	
	# output stats
	print getCurrentDateTime()."\t";
	print "Writing stat infos..\n";
	
	if ($sCount==1){$statOutFile=$sampleIDs[0]."_readStat.txt";}
	else {$statOutFile=$sCount."samples_readStat.txt";}
	
	open(TableOut , ">$statOutFile") || die "The file $statOutFile can not be openned.\n";
	print TableOut "SampleID\tCleanPE(#)\tJoined(#)\tQualified(#)\tNoChime(#)\tBase(nt)\tAvgLen(nt)\tAvgQ\tQ20%\tQ30%\tGC%\tEffective%\t";
	print TableOut "Mean\tMode\tMedian\tMin\tMax\tStdev\n";
	for ($i=0;$i< $sCount ;$i++)
	{
		print TableOut "$sampleIDs[$i]\t".commify($cleanReadCount{$sampleIDs[$i]})."\t".commify($joinedReadCount{$sampleIDs[$i]})."\t".commify($qualifiedReadCount{$sampleIDs[$i]})."\t".commify($effectiveReadCount{$sampleIDs[$i]})."\t".
		               commify($effectiveBaseCount{$sampleIDs[$i]})."\t".commify($effectiveAvgLen{$sampleIDs[$i]})."\t".commify($effectiveAvgQ{$sampleIDs[$i]})."\t".
		               "$effectiveQ20pct{$sampleIDs[$i]}\t$effectiveQ30pct{$sampleIDs[$i]}\t$effectiveGCpct{$sampleIDs[$i]}\t$effectivePct{$sampleIDs[$i]}\t";
	
		print TableOut commify($effectiveMeanLen{$sampleIDs[$i]})."\t".commify($effectiveModeLen{$sampleIDs[$i]})."\t".commify($effectiveMedianLen{$sampleIDs[$i]})."\t".commify($effectiveMinLen{$sampleIDs[$i]})."\t".commify($effectiveMaxLen{$sampleIDs[$i]})."\t".commify($effectiveStdevLen{$sampleIDs[$i]});
		print TableOut "\n";
				
	}
	close(TableOut);

}

print getCurrentDateTime()."\t";
print "Processing 16S analysis...\n";
# make merge groups
#$mergeFastaFile="mergegroups";
#if (-e "$folder[5]/$mergeFastaFile.fasta"){unlink("$folder[5]/$mergeFastaFile.fasta");}
chdir ($folder[5]);
print "\ncd $folder[5]\n";

for ($i=0;$i< $sCount ;$i++)
{
	print ($i+1);
	print ". $sampleIDs[$i]:\n";
	
	if ( $baseSpaceTag eq "Y" && !-e "../$folder[8]/$sampleIDs[$i]_S1_L001_R1_001.fastq.gz" ) #for illumina basespace 16S analysis app
	{
		print getCurrentDateTime()."\t";
		print "fastaQual2fastq.pl ../$folder[4]/$sampleIDs[$i].effective.fasta ../$folder[4]/$sampleIDs[$i].effective.qual 33 1000 &\n";
		`/export/EC1680U/software/anaconda2/bin/mothur "#make.fastq(fasta=../$folder[4]/$sampleIDs[$i].effective.fasta, qfile=../$folder[4]/$sampleIDs[$i].effective.qual)"`;
#		`/export/EC1680U/software/fastaQual2fastq.pl ../$folder[4]/$sampleIDs[$i].effective.fasta ../$folder[4]/$sampleIDs[$i].effective.qual 33 1000`;
		`gzip -c ../$folder[4]/$sampleIDs[$i].effecitive.fastq > ../$folder[8]/$sampleIDs[$i]_S1_L001_R1_001.fastq.gz`;
	}
	
	if (!-d "$sampleIDs[$i].$folderSuffix[5]") {mkdir "$sampleIDs[$i].$folderSuffix[5]";}

#	if (-e "$sampleIDs[$i].$folderSuffix[5]/$sampleIDs[$i].unique.filter.opti_mcc.0.03.biom"){ push @biomNames, "$sampleIDs[$i].biom"; next;}
	if(-e "../$folder[6]/$sampleIDs[$i].taxID.filter.dense.dedup.biom"){next;}
	
	chdir ("$sampleIDs[$i].$folderSuffix[5]");
	print "cd $sampleIDs[$i].$folderSuffix[5]\n";
	
	if (!-l "$sampleIDs[$i].$fastaSuffix[5]" && -e "../../$folder[4]/$sampleIDs[$i].$fastaSuffix[4]"){`ln -s ../../$folder[4]/$sampleIDs[$i].$fastaSuffix[4] $sampleIDs[$i].$fastaSuffix[5]`;}
	#if (!-l "$sampleIDs[$i].$qualSuffix[5]" && -e "../../$folder[4]/$sampleIDs[$i].$qualSuffix[4]"){`ln -s ../../$folder[4]/$sampleIDs[$i].$qualSuffix[4] $sampleIDs[$i].$qualSuffix[5]`;}
	
	print getCurrentDateTime()."\t";
	print "mothur \"#make.group(fasta=$sampleIDDash[$i].$fastaSuffix[5], groups=$sampleIDDash[$i])\"\n";
	if(! -e "$sampleIDs[$i].groups"){
	`/export/EC1680U/software/anaconda2/bin/mothur "#make.group(fasta=$sampleIDDash[$i].$fastaSuffix[5], groups=$sampleIDDash[$i])"`;
	}
	#Output File Names: 
	#T635-1TEST.groups
	print getCurrentDateTime()."\t";
	print "mothur \"#unique.seqs(fasta=$sampleIDs[$i].$fastaSuffix[5])\"\n";
	if(! -e "$sampleIDs[$i].unique.fasta"){
	`/export/EC1680U/software/anaconda2/bin/mothur "#unique.seqs(fasta=$sampleIDs[$i].$fastaSuffix[5])"`;
	}
	#Output File Names: 
	#T635-1TEST.names
	#T635-1TEST.unique.fasta
	print getCurrentDateTime()."\t";
	print "mothur \"#count.seqs(processors=4, name=$sampleIDs[$i].names, group=$sampleIDs[$i].groups)\"\n";
	if(! -e "$sampleIDs[$i].count_table"){
	`/export/EC1680U/software/anaconda2/bin/mothur "#count.seqs(processors=4, name=$sampleIDs[$i].names, group=$sampleIDs[$i].groups)"`;
	}
	#Output File Names: 
	#T635-1TEST.count_table
	
	if (!-e "$sampleIDs[$i].unique.raw_taxonomy.wang.taxonomy")
	{
	print getCurrentDateTime()."\t";
	print "mothur \"#classify.seqs(processors=30,fasta=$sampleIDDash[$i].unique.fasta,name=$sampleIDDash[$i].names,reference=$refFastaDB[1],taxonomy=$refTaxonomy[1], cutoff=0, method=wang, probs=f)\"\n";
	`/export/EC1680U/software/anaconda2/bin/mothur "#classify.seqs(processors=30,fasta=$sampleIDDash[$i].unique.fasta,name=$sampleIDDash[$i].names,reference=$refFastaDB[1],taxonomy=$refTaxonomy[1], cutoff=0, method=wang, probs=f)"`;
	
	#It took 1251 secs to classify 127500 sequences.
	#
	#Reading T635-1TEST.names...  Done.
	#
	#It took 4 secs to create the summary file for 127500 sequences.
	#
	#
	#Output File Names: 
	#T635-1TEST.unique.raw_taxonomy.wang.taxonomy
	#T635-1TEST.unique.raw_taxonomy.wang.tax.summary
	#T635-1TEST.unique.raw_taxonomy.wang.flip.accnos
	}


$removeLineageTag=0;
if ($removeLineageTag==1) 
{
	print getCurrentDateTime()."\t";
	print "mothur \"#remove.lineage(taxonomy=$sampleIDs[$i].unique.raw_taxonomy.wang.taxonomy,taxon=unknown,fasta=$sampleIDs[$i].unique.fasta,name=$sampleIDs[$i].names, group=$sampleIDs[$i].groups)\"\n";
	`/export/EC1680U/software/anaconda2/bin/mothur "#remove.lineage(taxonomy=$sampleIDs[$i].unique.raw_taxonomy.wang.taxonomy,taxon=Chloroplast-Mitochondria-unknown-mitochondria-chloroplast,fasta=$sampleIDs[$i].unique.fasta,name=$sampleIDs[$i].names, group=$sampleIDs[$i].groups)"`;

	#Output File Names: 
	#T635-1TEST.unique.raw_taxonomy.wang.pick.taxonomy
	#T635-1TEST.pick.names
	#T635-1TEST.unique.pick.fasta
	#T635-1TEST.pick.groups
}

	if (!-e "$sampleIDs[$i].unique.align")
	{
	print getCurrentDateTime()."\t";
	print "mothur \"#align.seqs(processors=30, flip=F, save=T, candidate=$sampleIDDash[$i].unique.fasta, template=$refAlignDB[1])\"\n";
	`/export/EC1680U/software/anaconda2/bin/mothur "#align.seqs(processors=30, flip=F, candidate=$sampleIDDash[$i].unique.fasta, template=$refAlignDB[1])"`;

	#	It took 1130 secs to align 127500 sequences.
	#Output File Names: 
	#T635-1TEST.unique.align
	#T635-1TEST.unique.align.report
	#T635-1TEST.unique.flip.accnos
	}
	
	print getCurrentDateTime()."\t";
	if(! -e "$sampleIDs[$i].unique.summary"){
	print "mothur \"#summary.seqs(processors=4,fasta=$sampleIDs[$i].unique.align, count=$sampleIDs[$i].count_table)\"\n";
	`/export/EC1680U/software/anaconda2/bin/mothur "#summary.seqs(processors=4,fasta=$sampleIDs[$i].unique.align, count=$sampleIDs[$i].count_table)"`;
	}
	#Output File Names: 
	#T635-1TEST.unique.summary

	print getCurrentDateTime()."\t";
	if(! -e "$sampleIDs[$i].unique.filter.fasta"){
	print "mothur \"#filter.seqs(processors=10,fasta=$sampleIDDash[$i].unique.align, vertical=T)\"\n";
	`/export/EC1680U/software/anaconda2/bin/mothur "#filter.seqs(processors=10,fasta=$sampleIDDash[$i].unique.align, vertical=T)"`;
	}
	#Length of filtered alignment: 1306
	#Number of columns removed: 6376
	#Length of the original alignment: 7682
	#Number of sequences used to construct filter: 1713212
	#
	#Output File Names: 
	#T635-1TEST.filter
	#T635-1TEST.unique.filter.fasta

	if (!-e "$sampleIDs[$i].unique.filter.dist")
	{
	#dist need same length alignment 
	print getCurrentDateTime()."\t";   
	print "mothur \"#dist.seqs(processors=30,fasta=$sampleIDs[$i].unique.filter.fasta,calc=onegap,countends=F,cutoff=0.03)\"\n";
	`/export/EC1680U/software/anaconda2/bin/mothur "#dist.seqs(processors=30,fasta=$sampleIDs[$i].unique.filter.fasta,calc=onegap,countends=F,cutoff=0.03)"`;
	#Output File Names: 
	#T635-1TEST.unique.filter.dist
	#
	#It took 2105 to calculate the distances for 127500 sequences.
	}
	if (!-e "$sampleIDs[$i].unique.filter.opti_mcc.list")
	{
	print getCurrentDateTime()."\t";
	print "mothur \"#cluster(column=$sampleIDs[$i].unique.filter.dist, name=$sampleIDs[$i].names, cutoff=0.03, method=opti)\"\n";
	`/export/EC1680U/software/anaconda2/bin/mothur "#cluster(column=$sampleIDs[$i].unique.filter.dist, name=$sampleIDs[$i].names, cutoff=0.03, method=opti)"`; #By default cluster() uses the average neighbor algorithm;

	#changed cutoff to 0
	#
	#Output File Names: 
	#T635-1TEST.unique.filter.opti_mcc.sabund
	#T635-1TEST.unique.filter.opti_mcc.rabund
	#T635-1TEST.unique.filter.opti_mcc.list
	#
	#It took 36 seconds to cluster
	}	
	print getCurrentDateTime()."\t";
	if(! -e "$sampleIDs[$i].unique.filter.opti_mcc.shared"){
	print "mothur \"#make.shared(list=$sampleIDs[$i].unique.filter.opti_mcc.list, group=$sampleIDs[$i].groups)\"\n";
	`/export/EC1680U/software/anaconda2/bin/mothur "#make.shared(list=$sampleIDs[$i].unique.filter.opti_mcc.list, group=$sampleIDs[$i].groups)"`;
	}
	#Output File Names: 
	#T635-1TEST.unique.filter.opti_mcc.shared
	#T635-1TEST.unique.filter.opti_mcc.T635-1TEST.rabund

	print getCurrentDateTime()."\t";
	if(! -e "$sampleIDs[$i].unique.filter.opti_mcc.0.03.cons.taxonomy"){
	print "mothur \"#classify.otu(list=$sampleIDs[$i].unique.filter.opti_mcc.list,count=$sampleIDs[$i].count_table, taxonomy=$sampleIDs[$i].unique.raw_taxonomy.wang.taxonomy, cutoff=80, label=0.03, probs=f)\"\n";
	`/export/EC1680U/software/anaconda2/bin/mothur "#classify.otu(list=$sampleIDs[$i].unique.filter.opti_mcc.list,count=$sampleIDs[$i].count_table, taxonomy=$sampleIDs[$i].unique.raw_taxonomy.wang.taxonomy, cutoff=80, label=0.03, probs=f)"`;
	}

	#Output File Names: 
	#$sampleIDs[$i].unique.filter.opti_mcc.0.03.cons.taxonomy
	#$sampleIDs[$i].unique.filter.opti_mcc.0.03.cons.tax.summary

	print getCurrentDateTime()."\t";
	if(! -e "$sampleIDs[$i].unique.filter.opti_mcc.0.03.biom"){
	print "mothur \"#make.biom(shared=$sampleIDs[$i].unique.filter.opti_mcc.shared,constaxonomy=$sampleIDs[$i].unique.filter.opti_mcc.0.03.cons.taxonomy, matrixtype=dense, reftaxonomy=$refTaxonomy[1]\"\n";
	`/export/EC1680U/software/anaconda2/bin/mothur "#make.biom(shared=$sampleIDs[$i].unique.filter.opti_mcc.shared,constaxonomy=$sampleIDs[$i].unique.filter.opti_mcc.0.03.cons.taxonomy, matrixtype=dense, reftaxonomy=$refTaxonomy[1])"`;
	}
	#Output File Names: 
	#$sampleIDs[$i].unique.filter.opti_mcc.0.03.biom
	
	print getCurrentDateTime()."\t";
	if(! -e "$sampleIDs[$i].unique.filter.opti_mcc.0.03.rep.fasta"){
	print "mothur \"#get.oturep(column=$sampleIDs[$i].unique.filter.dist, list=$sampleIDs[$i].unique.filter.opti_mcc.list, name=$sampleIDs[$i].names, fasta=$sampleIDs[$i].unique.fasta)\"\n";
	`/export/EC1680U/software/anaconda2/bin/mothur "#get.oturep(column=$sampleIDs[$i].unique.filter.dist, list=$sampleIDs[$i].unique.filter.opti_mcc.list, name=$sampleIDs[$i].names, fasta=$sampleIDs[$i].unique.fasta)"`;
	}
	#Output File Names: 
	#C1315A.unique.filter.opti_mcc.0.03.rep.names
	#C1315A.unique.filter.opti_mcc.0.03.rep.fasta

	#print getCurrentDateTime()."\t";
	#print "mothur \"#tree.shared(processors=12, shared=$sampleIDs[$i].unique.filter.opti_mcc.shared, name=$sampleIDs[$i].names, groups=all)\"\n";
	#`mothur "#tree.shared(processors=12, shared=$sampleIDs[$i].unique.filter.opti_mcc.shared, name=$sampleIDs[$i].names, groups=all)"`;
	#Output File Names: 
	#You have not provided enough valid groups.  I cannot run the command. ==> less than 2 groups
	
	
	chdir ("../../$folder[6]");
	print "cd ../../$folder[6]\n";
	`pwd`;
	if (! -e "$sampleIDs[$i].otuID.raw.biom")
	{
		print "cp ../$folder[5]/$sampleIDs[$i].$folderSuffix[5]/$sampleIDs[$i].unique.filter.opti_mcc.0.03.biom $sampleIDs[$i].otuID.raw.biom\n";
		`cp ../$folder[5]/$sampleIDs[$i].$folderSuffix[5]/$sampleIDs[$i].unique.filter.opti_mcc.0.03.biom $sampleIDs[$i].otuID.raw.biom`;
	}


	##### Merge otuID biom files and pipe with JSYF pipeline #####
	if(! -e "$sampleIDs[$i].otuID.filter.biom"){
	system("$qiime_docker filter_otus_from_otu_table.py -i $sampleIDs[$i].otuID.raw.biom -o $sampleIDs[$i].otuID.filter.biom --min_count_fraction 0.00005"); # ref: https://goo.gl/EMwR9N
	}

	## Convert otuID.filter to taxID.filter
	if(! -e "$sampleIDs[$i].results.blast6"){
	system("/export/EC1680U/software/anaconda2/bin/vsearch ".
	    "--usearch_global ../$folder[5]/$sampleIDs[$i].$folderSuffix[5]/$sampleIDs[$i].unique.filter.opti_mcc.0.03.rep.fasta ".
		"--db $refFastaDB[1] --id 0.1 --iddef 1 --blast6out $sampleIDs[$i].results.blast6 --maxaccepts 1 --threads 32");
	}

	## Powered by JSYF
	system("perl /export/EC1680U/perl/bin/16S/biomProc.pl sparseToDense ".
	    "$sampleIDs[$i].otuID.filter.biom ".
	    "$sampleIDs[$i].otuID.filter.dense.biom");

	system("perl /export/EC1680U/perl/bin/16S/biomProc.pl idChangeUsingFasta ".
	    "$sampleIDs[$i].otuID.filter.dense.biom ".
    	"$sampleIDs[$i].taxID.filter.dense.biom ".
	    "../$folder[5]/$sampleIDs[$i].$folderSuffix[5]/$sampleIDs[$i].unique.filter.opti_mcc.0.03.rep.fasta,$sampleIDs[$i].results.blast6");

	system("perl /export/EC1680U/perl/bin/16S/biomProc.pl uniqById ".
	    "$sampleIDs[$i].taxID.filter.dense.biom ".
    	"$sampleIDs[$i].taxID.filter.dense.dedup.biom");
	###
	
	system("ln -s $sampleIDs[$i].taxID.filter.dense.dedup.biom $sampleIDs[$i].biom");
	
	##############################################################
	
	chdir ("../$folder[5]");
	print "cd ../$folder[5]\n";
	push @biomNames, "$sampleIDs[$i].biom";
	#push @groupNames, $sampleIDDash[$i];
	
	#`cat $folder[4]/$sampleIDs[$i].$fastaSuffix[4] >> $folder[5]/$mergeFastaFile.fasta`;
	#$fastaNames= join('-',@fastaNames);
	#$groupNames= join('-',@groupNames);

	#chdir ("..");
}
chdir ("..");
print "cd ..\n";

##### Merge all biom files #####
chdir("6_biom");
my $all_biom = `ls *.taxID.filter.dense.dedup.biom | sed ':a;N;\$!ba;s/\\n/,/g'`; chomp($all_biom);
system("$qiime_docker merge_otu_tables.py -i $all_biom -o otu_table.biom");
chdir("..");
#####

#$OTUdataFolder = $folder[7]/OTUs_eachSample = 16S_OUTPUT/OTUs_eachSample
if ($OTUstatTag eq "Y")
{
	if (!-d $OTUdataFolder){mkdir $OTUdataFolder;}
	chdir ($OTUdataFolder);
	
	if ($sCount==1){$statOTUfile=$sampleIDs[0].".OTUstat.txt";}
	else {$statOTUfile=$sCount."samples.OTUstat.txt";}
	if ($sCount==1){$statTaxaFile=$sampleIDs[0].".TAXONOMYstat.txt";}
	else {$statTaxaFile=$sCount."samples.TAXONOMYstat.txt";}
	
	
	
	
	print "\n\nStat. the effective reads, bp, Taxon Reads, Unique Reads, OTU number of each sample...\n\n";
	print getCurrentDateTime()."\t";
	print "Analysing OTU infos..\n";
	# count reads and stats
	for ($i=0;$i< $sCount ;$i++)
	{
		if (!-d "$sampleIDs[$i]"){mkdir ($sampleIDs[$i]);}
		chdir ($sampleIDs[$i]);
		if (!-l "$sampleIDs[$i].rep_set.fasta")
		{
			print "ln -s ../../../$folder[5]/$sampleIDs[$i].$folderSuffix[5]/$sampleIDs[$i].unique.filter.opti_mcc.0.03.rep.fasta $sampleIDs[$i].rep_set.fasta\n";
			`ln -s ../../../$folder[5]/$sampleIDs[$i].$folderSuffix[5]/$sampleIDs[$i].unique.filter.opti_mcc.0.03.rep.fasta $sampleIDs[$i].rep_set.fasta`;
		}
		if (!-l "$sampleIDs[$i].rep_set.names")
		{
			print "ln -s ../../../$folder[5]/$sampleIDs[$i].$folderSuffix[5]/$sampleIDs[$i].unique.filter.opti_mcc.0.03.rep.names $sampleIDs[$i].rep_set.names\n";
			`ln -s ../../../$folder[5]/$sampleIDs[$i].$folderSuffix[5]/$sampleIDs[$i].unique.filter.opti_mcc.0.03.rep.names $sampleIDs[$i].rep_set.names`;
		}
		if (!-l "$sampleIDs[$i].rep_set.taxonomy")
		{
			print "ln -s ../../../$folder[5]/$sampleIDs[$i].$folderSuffix[5]/$sampleIDs[$i].unique.filter.opti_mcc.0.03.cons.taxonomy $sampleIDs[$i].rep_set.taxonomy\n";
			`ln -s ../../../$folder[5]/$sampleIDs[$i].$folderSuffix[5]/$sampleIDs[$i].unique.filter.opti_mcc.0.03.cons.taxonomy $sampleIDs[$i].rep_set.taxonomy`;
		}
		
		### 2017/11/14 added by kentchen
		`/export/EC1680U/perl/bin/tab2xlsx.pl ../../../$folder[5]/$sampleIDs[$i].$folderSuffix[5]/$sampleIDs[$i].unique.align.report $sampleIDs[$i].unique.align.report.xlsx`;

		if (!defined($effectiveBaseCount{$sampleIDs[$i]}) || $effectiveBaseCount{$sampleIDs[$i]}==0)
		{
				if (-e "../../../$folder[4]/$sampleIDs[$i].$fastaSuffix[4].gz")
				{($total,$totalLen,$avgLen,$GCpct,$mean,$mode,$median,$min,$max,$stdev)=statFasta("../../../$folder[4]/$sampleIDs[$i].$fastaSuffix[4].gz",1);}
				elsif (-e "../../../$folder[4]/$sampleIDs[$i].$fastaSuffix[4]")
				{($total,$totalLen,$avgLen,$GCpct,$mean,$mode,$median,$min,$max,$stdev)=statFasta("../../../$folder[4]/$sampleIDs[$i].$fastaSuffix[4]");}
			
				$effectiveReadCount{$sampleIDs[$i]}=$total;
				$effectiveBaseCount{$sampleIDs[$i]}=$totalLen;

		}

		($taxonReads,$NrReads,$OTUcount,@taxonComposition)=statOTUs($sampleIDs[$i],"../../../$folder[5]/$sampleIDs[$i].$folderSuffix[5]/$sampleIDs[$i].unique.fasta","$sampleIDs[$i].rep_set.taxonomy");
		#(@taxonComposition,@taxonPercentage)=statTAXONs($sampleIDs[$i],"$sampleIDs[$i].rep_set.taxonomy");
			
		
		$OTUstatStr[$i]="$sampleIDs[$i]\t".commify($effectiveReadCount{$sampleIDs[$i]})."\t".commify($effectiveBaseCount{$sampleIDs[$i]})."\t".commify($taxonReads)."\t".commify($NrReads)."\t".commify($OTUcount);
	

		@taxonPercentage = ( sprintf("%.1f%%", $taxonComposition[0]/$effectiveReadCount{$sampleIDs[$i]}*100), sprintf("%.1f%%", $taxonComposition[1]/$effectiveReadCount{$sampleIDs[$i]}*100),sprintf("%.1f%%", $taxonComposition[2]/$effectiveReadCount{$sampleIDs[$i]}*100),sprintf("%.1f%%", $taxonComposition[3]/$effectiveReadCount{$sampleIDs[$i]}*100),sprintf("%.1f%%", $taxonComposition[4]/$effectiveReadCount{$sampleIDs[$i]}*100),sprintf("%.1f%%", $taxonComposition[5]/$effectiveReadCount{$sampleIDs[$i]}*100),sprintf("%.1f%%", $taxonComposition[6]/$effectiveReadCount{$sampleIDs[$i]}*100) );
		
		
		$TAXAstatStr[$i]="$sampleIDs[$i]\tAmount\t".join("\t",(map {commify($_)} @taxonComposition))."\n".
										"\tPercentage(%)\t".join("\t",@taxonPercentage)."\n";
	
		chdir ("..");
	}
	

	
	open(TableOut , ">$statOTUfile") || die "The file $statOTUfile can not be openned.\n";
	print TableOut "SampleID\tEffectiveReads(#)\tEffectiveBases(bp)\tTaxonReads(#)\tNrReads(#)\tOTUs(#)\n";
	print TableOut join("\n",@OTUstatStr);
	print TableOut "\n";
	close TableOut;
	
	
	
	open(TableOut , ">$statTaxaFile") || die "The file $statTaxaFile can not be openned.\n";
	print TableOut "SampleID\tKindom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
	print TableOut join("",@TAXAstatStr);
	close TableOut;
	
	chdir ("../..");
}


	$biomOutFileSH="runMergeBioms.sh";
	open(shBiomOut , ">$folder[6]/$biomOutFileSH") || die "The file $folder[6]/$biomOutFileSH can not be openned.\n";
	$jointBiomNames=join ",",@biomNames;
	print shBiomOut "merge_otu_tables.py -i $jointBiomNames -o otu_table.biom\n";
#View statistics of the OTU table
	print shBiomOut "statistics of the OTU table\n";
	print shBiomOut "biom summarize-table -i otu_table.biom -o otu_table_summary.txt\n";

#Make Taxa Summary Charts
	print shBiomOut "echo 'Summarize taxa'\n";
	print shBiomOut "rm -rf taxa_summary\n";
# don't use mapFile here for summarize_taxa
	print shBiomOut "summarize_taxa.py -i otu_table.biom -o taxa_plots -a  -L 2,3,4,5,6,7,8\n";
	print shBiomOut "cd taxa_plots\n";
	print shBiomOut "plot_taxa_summary.py -i otu_table_L2.txt,otu_table_L3.txt,otu_table_L4.txt,otu_table_L5.txt,otu_table_L6.txt,otu_table_L7.txt,otu_table_L8.txt -l Kindom,Phylum,Class,Order,Family,Genus,Species -o .\n";
	print shBiomOut "mv -f otu* raw_data/. \n";
	print shBiomOut "cd ..\n";
	print shBiomOut "mv taxa_plots ../$folder[7]\n";
#OTU Heatmap
	print shBiomOut "echo 'OTU Heatmap'\n";
	print shBiomOut "make_otu_heatmap.py -i otu_table.biom -o OTUs_heatmap\n";
	print shBiomOut "make_otu_heatmap_html.py -i otu_table.biom -o OTUs_heatmap/  \n";
	print shBiomOut "mv OTUs_heatmap ../$folder[7]\n";

#plot rank-abundance curve
	print shBiomOut "echo 'Rank-abundance curve'\n";
	print shBiomOut "plot_rank_abundance_graph.py -i otu_table.biom -s '*' -o rank_abundance_graph.L -x -f pdf\n";
	print shBiomOut "plot_rank_abundance_graph.py -i otu_table.biom -s '*' -o rank_abundance_graph.L -x -f png\n";
	print shBiomOut "plot_rank_abundance_graph.py -i otu_table.biom -s '*' -o rank_abundance_graph.L -x -f svg \n";
# no legends
	print shBiomOut "plot_rank_abundance_graph.py -i otu_table.biom -s '*' -o rank_abundance_graph -n -x -f pdf\n";
	print shBiomOut "plot_rank_abundance_graph.py -i otu_table.biom -s '*' -o rank_abundance_graph -n -x -f png\n";
	print shBiomOut "plot_rank_abundance_graph.py -i otu_table.biom -s '*' -o rank_abundance_graph -n -x -f svg \n";
	print shBiomOut "../$folder[7]/rank_abundance_plots\n";
	print shBiomOut "mv rank_abundance_graph* ../$folder[7]/rank_abundance_plots/.\n";
	print shBiomOut "rm -rf ../$folder[7]/rank_abundance_plots/*.txt\n";
	
	
#Calculate alpha diversity on each sample in an otu table, using a variety of alpha diversity metrics
#Known metrics are: ACE, berger_parker_d, brillouin_d, chao1, chao1_confidence, dominance, doubles, enspie, equitability, esty_ci, fisher_alpha, gini_index, 
#goods_coverage, heip_e, kempton_taylor_q, margalef, mcintosh_d, mcintosh_e, menhinick, michaelis_menten_fit, observed_species, osd, simpson_reciprocal, 
#robbins, shannon, simpson, simpson_e, singles, strong, PD_whole_tree

	print shBiomOut "echo 'Alpha diversity'\n";
	print shBiomOut "echo 'alpha_diversity:metrics chao1,ACE,shannon,simpson_reciprocal,fisher_alpha,goods_coverage,observed_species' > alpha_params.txt\n";
	print shBiomOut "alpha_diversity.py -i otu_tables/ -o alpha_div/ -m chao1,ACE,shannon,simpson_reciprocal,fisher_alpha,goods_coverage,observed_species\n";
	print shBiomOut "collate_alpha.py -i alpha_div/ -o alpha_div_collated/\n";
	print shBiomOut "make_rarefaction_plots.py -i alpha_div_collated/ -o rarefaction_plots -g pdf --generate_average_tables\n";
	print shBiomOut "alpha_rarefaction.py  -i otu_table.biom -p alpha_params.txt -m Fasting_Map.txt -o alpha_rarefaction -p alpha_params.txt -t otus/rep_set.tre \n";



#OTU Network
	print shBiomOut "echo 'OTU Network'\n";
	print shBiomOut "make_otu_network.py -m mapFile.txt -i otu_table.biom -o OTU_Network\n";

	print shBiomOut "tree -F -f --nolinks --dirsfirst --nolinks -H  meta16S_OUTPUT meta16S_OUTPUT -o meta16S_OUTPUT/dirTree.html\n";


	close(shBiomOut);
	
	
print getCurrentDateTime()."\t";
print "Done $0 !!!\n\n";
print "Check the following result files:\n";
print "$OTUdataFolder/$statOTUfile\n";
print "$OTUdataFolder/$statTaxaFile\n";

print "Please run the script: '$folder[6]/$biomOutFileSH' under QIIME (Unbuntu VirtualBox)\n\n";

exit();

#####Subrotines######

sub decimalcommas  
{
	local $_  = shift;
	1 while s/^([-+]?\d+)(\d{3})/$1,$2/;
	return $_;
}
sub cmplc # complement a string using IUPAC codes for [deoxy]ribonucleotides
{
	local($_) = shift(@_);
	$iupac_nt   = "ACGTUMRWSYKVHDBNacgtumrwsykvhdbn[]" ;
	$iupac_cmpl = "TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn][" ;
	eval "tr/53 $iupac_nt/35 $iupac_cmpl/";
	return $_ ;
}

sub cmplr # simply reverse a string character-for-character
{
	local($seq) = shift(@_);
	local($rev) = scalar(reverse($seq));
	return $rev ;
}




sub trim
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub commify
{
	my $number = shift;
	$number =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
	return $number;
}


sub timeStamp
{
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime(); 
	my $yr = $year+1900;
	my $month = $mon+1;
	$month = "0$month" if $month<10;
	$day = "0$day" if $day<10;
	$hour = "0$hour" if $hour<10;
	$min = "0$min" if $min<10;
	$sec = "0$sec" if $sec<10;

	$timeDec =  "$yr$month$day$hour$min$sec";
	return $timeDec;
}


sub statFastq
{
	my ($fastq,$gzTag) = @_;
	my $statStr;
	
	#print "$fastq,$gzTag\n";
	if (!defined($gzTag) || $gzTag==0)
	{
		$statStr=`cat $fastq | awk '((NR-2)%4==0){ read=\$1;rlen=length(read);total++;totalLen+=rlen }END{ print total,totalLen,totalLen/total }'`;
	}
	else
	{
		$statStr=`gzip -dc $fastq | awk '((NR-2)%4==0){ read=\$1;rlen=length(read);total++;totalLen+=rlen }END{ print total,totalLen,totalLen/total }'`;
		#print "gzip -dc $fastq | awk '((NR-2)%4==0){ read=\$1;rlen=length(read);total++;totalLen+=rlen }END{ print total,totalLen,totalLen/total }'\n";
	}
	
	chomp $statStr;
	my @items=split / +/,$statStr;
	
#`cat $fastq | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};
#print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'`

	return @items;
	#return ($items[0],$items[1],$items[2]); #($total,$totalLen,$avgLen)
}

sub statFasta
{
	#my $fasta = shift;
	my ($fasta,$gzTag) = @_;
	my ($seqName,%seqLength,%seqSequence);
	my ($total,$totalLen,$avgLen,$GCpct);
	my (@lines,$line,$GC);
	
	if (!defined($gzTag) || $gzTag==0)
	{
		@lines=`cat $fasta`;
	}
	else
	{
		@lines=`gzip -dc $fasta`;
	}

	#open (file_in,$fasta)||die "cannot open this file $fasta";
	#while(defined($line=<file_in>)) 
	foreach $line (@lines)
	{
		chomp $line;
		#$line =~ s/\r\n?//g;
		$line=trim($line);
		
		if ($line eq "" || $line =~ /^ +$/){next;}
		if ($line =~ /^\>(.+)$/)
		{ 
			$seqName=$1;
			$seqLength{$seqName}=0;
			$seqSequence{$seqName}="";
			$total++;
		}
		else
		{
			$seqSequence{$seqName}.=$line;
		  $seqLength{$seqName}+=length($line);
		  $totalLen+=length($line);
		  $GC++ while $line =~ /[GC]/ig;
		}
	}
	#close(file_in);

# calculate median


my @seqLengthTopDown=  sort { $seqLength{$b} <=> $seqLength{$a} } keys %seqLength;

my $piddle = pdl @seqLengthTopDown;
my ($mean,$prms,$median,$min,$max,$adev,$stdev) = statsover $piddle;
my $mode = mode($piddle);

	
	$avgLen=sprintf("%.02f", $totalLen/$total);
	$GCpct=sprintf("%.02f", $GC/$totalLen*100);
	
	return ($total,$totalLen,$avgLen,$GCpct,$mean,$mode,$median,$min,$max,$stdev);
}

sub statQual
{
	#my $qual = shift;
	my ($qual,$gzTag) = @_;
	my ($Q20pct,$Q30pct,$avgQ);
	my (@lines,$line,$ctgName,@pQarray,$pQ,$pQ20,$pQ30,$pQtotal,$pQsum);
	$pQ=$pQ20=$pQ30=0;
	$pQtotal=$pQsum=0;
	
	if (!defined($gzTag) || $gzTag==0)
	{
		@lines=`cat $qual`;
	}
	else
	{
		@lines=`gzip -dc $qual`;
	}
	
	#open (file_in,$qual)||die "cannot open this file $qual";
	#while(defined($line=<file_in>)) 
	foreach $line (@lines)
	{
		chomp $line;
		$line =~ s/\r\n?//g;
		$line=trim($line);
		
		if ($line eq "" || $line =~ /^ +$/){next;}

		if ($line =~ /^(.)?\>(.+)$/)
		{ 
			$ctgName=$1;
		}
		else
		{
		  @pQarray=split / +/,$line;
			
			#print "pQtotal|pQsum|pQ20|pQ30=$pQtotal|$pQsum|$pQ20|$pQ30\n";
			foreach $pQ (@pQarray)
			{
				if ($pQ !~ /^[0-9]+$/){next;} # in case of strange charactors would exist in the end of the qual file from the mothur instruction: fastq.info
				
				$pQtotal++;
				$pQsum+=$pQ; #
				if ($pQ >= 20){$pQ20++;}
				if ($pQ >= 30){$pQ30++;}
				
			}
			
		}
	}
	#close(file_in);
	
	$avgQ= sprintf("%.02f",$pQsum/$pQtotal);
	$Q20pct= sprintf("%.02f", $pQ20/$pQtotal*100);
	$Q30pct= sprintf("%.02f", $pQ30/$pQtotal*100);
	
	#print "($pQtotal,$avgQ, $Q20pct, $Q30pct)\n";
	return ($avgQ, $Q20pct, $Q30pct);
}

sub getCurrentDateTime 
{

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $nice_timestamp = sprintf ( "%04d-%02d-%02d %02d:%02d:%02d",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return "[$nice_timestamp]";
}



#($taxonReads,$NrReads,$OTUcount)=statOTUs($sampleIDs[$i],"$sampleIDs[$i].rep_set.fasta","../../$folder[5]/$sampleIDs[$i].$folderSuffix[5]/$sampleIDs[$i].unique.fasta","$sampleIDs[$i].rep_set.taxonomy");

#>M02662_7_000000000-A9RN7_1_1101_8319_7543	Otu00001|342

sub statOTUs
{

	my ($sampleID,$UNQfasta,$taxonomy) = @_;
	

	if (!-e $UNQfasta){die "$UNQfasta not found!!\n\n"; }
	if (!-e $taxonomy){die "$taxonomy not found!!\n\n"; }
	
	my ($taxonReads,$NrReads,$OTUcount)=(0,0,0);
	my (@lines,$line);
	my (%OTUrepReadNum,$OTUtaxonStr,@OTUtaxons,$OTUid);
	my ($nK,$nP,$nC,$nO,$nF,$nG,$nS)=(0,0,0,0,0,0,0);
	
	$NrReads=`cat $UNQfasta | grep -c '>'`; chomp $NrReads;
	
	@lines=`cat $taxonomy`;
	foreach $line (@lines)
	{
		chomp $line;
		
		if ($line =~ /^OTU.+Taxonomy$/){next;}
		
		if ($line =~ /^(Otu[0-9]+)\t([0-9]+)\t(.+)$/)
		{
			$OTUid=$1;
			$OTUrepReadNum{$OTUid}=$2;
			$OTUtaxonStr=$3;  if (!defined($OTUtaxonStr) || $OTUtaxonStr eq ""){die "Cannot parse $taxonomy: '$line'\n\n";}
			$OTUcount++;
			
			if ($OTUtaxonStr !~/^Root;unclassified/){$taxonReads+=$OTUrepReadNum{$OTUid};}
			
			@OTUtaxons=split /;/,$OTUtaxonStr;
			#if ($OTUtaxons[0] =~/^Root/){}
			if (defined($OTUtaxons[0]) && $OTUtaxons[0] !~/unclassified/){$nK+=$OTUrepReadNum{$OTUid};}
			if (defined($OTUtaxons[1]) && $OTUtaxons[1] !~/unclassified/){$nP+=$OTUrepReadNum{$OTUid};}
			if (defined($OTUtaxons[2]) && $OTUtaxons[2] !~/unclassified/){$nC+=$OTUrepReadNum{$OTUid};}
			if (defined($OTUtaxons[3]) && $OTUtaxons[3] !~/unclassified/){$nO+=$OTUrepReadNum{$OTUid};}
			if (defined($OTUtaxons[4]) && $OTUtaxons[4] !~/unclassified/){$nF+=$OTUrepReadNum{$OTUid};}
			if (defined($OTUtaxons[5]) && $OTUtaxons[5] !~/unclassified/){$nG+=$OTUrepReadNum{$OTUid};}
			if (defined($OTUtaxons[6]) && $OTUtaxons[6] !~/unclassified/){$nS+=$OTUrepReadNum{$OTUid};}
			
		}
	}
	#@taxonComposition = (commify($nK),commify($nP),commify($nC),commify($nO),commify($nF),commify($nG),commify($nS));
	@taxonComposition = ($nK, $nP, $nC, $nO, $nF, $nG, $nS);
	
	return ($taxonReads,$NrReads,$OTUcount,@taxonComposition,);
}

