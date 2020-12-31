#! /usr/bin/env perl

# Written by Keh-Ming Wu, Ver 1.0, 2014.06.26
#
# 1. Dual-barcode demultiplexing for PE fastq reads (split fastq and remove barcode sequences from R1&R2 reads)
$version="Ver. 1.2";
$release_date="2014-08-11";
# 2. add single-barcode demultiplexing for PE fastq reads (split fastq and remove barcode sequences from R1 reads)
$version="Ver. 1.4";
$release_date="2014-09-28";

$programer = "Keh-Ming(Kevin) Wu";

# barcodes should be the same length in one barcode-infoSheet, or you should seperate them into different sheets and run for each sheet!~
# dualbarcodes	ATTACTCG	TCTCGCGC	C804-3A
# dualbarcodes	TCCGGAGA	TCTCGCGC	C804-3B             
# barcode	AAGTCTTG	PCR1
# barcode	AAGTTCCA	PCR2



use String::Approx 'amatch';
# use String::Approx 'aindex';
# use String::Approx 'arindex';
# use String::Approx 'aslice';
 
 
#use Getopt::Std;
#use strict;
use FileHandle;
STDOUT->autoflush(1);

##===================================================
MainPrg:

$Usage =	
    "=======================================\n".
		"* Something Wrong with the parameters *\n".
		"=======================================\n\n".
		"Usage:\n".
		"$0 R1.fq R2.fq barcodeFile outputFoler\n".
		"$0 140623-96xAmplicon_S1_L001_R1_001.clean.fastq 140623-96xAmplicon_S1_L001_R2_001.clean.fastq MS140017_barcodes.txt _demultiplexed _baseSpace\n\n";

if ($#ARGV < 1) { die "$0 $version ($release_date) by $programer\n$Usage\n";} 



if (defined($ARGV[0]))
{
	$fq[1]=$ARGV[0];
	if (!-e $fq[1]){die "$fq[1] not found!!\n\n";} else {$Rnum=1;}

	if ($fq[1]=~ /(.(clean|raw)\.(fastq|fq))/i){ $fqsuffix=$1;} 
	else 
	{
		$fqsuffix=".fastq";
	}
}
if (defined($ARGV[1]))
{
	$fq[2]=$ARGV[1];
	if (!-e $fq[2]){die "$fq[2] not found!!\n\n";} else {$Rnum=2;}
}
if (defined($ARGV[2]))
{
	$bFile=$ARGV[2];
	if (!-e $bFile){die "$bFile not found!!\n\n";}
}
else { die "barcode file not defined!!\n\n"; }

if (defined($ARGV[3])) {$outputFoler=$ARGV[3];}else{$outputFoler="_demultiplexed";}
if (!-d $outputFoler){`mkdir $outputFoler`;}

if (defined($ARGV[4])){$outputBSFoler=$ARGV[4]; if (!-d $outputBSFoler){`mkdir $outputBSFoler`;}}#else{$outputBSFoler="0_cleanFastqGz2BaseSpace";}



$date1=`date`;

$indexStatFile ="indexStat.$bFile";
open indexStat, ">$outputFoler/$indexStatFile" or die $!;
print indexStat $date1;
print indexStat "$0 ".join(" ", @ARGV)."\n\n";

#====================
$sampleIDerr=0;$barcodeErr=0;$bcCount=0;
$barcodeDupCheck=0;

print "Reading barcode File.... $bFile\n";
open(IN, $bFile) || die "The file $bFile can not be openned.\n";
while(defined($line=<IN>))
{ 
	chomp $line;  $line=~ s/\r\n?/\n/g;
	next if ($line =~ /^\#/);
	next if ($line eq "");
	
	@bcArray=split /\t/,$line;
	$itemNum=$#bcArray +1; 
	
	for ($i=0;$i< $itemNum ;$i++)
	{
		$bcArray[$i]=trim($bcArray[$i]);
	}
	
	$bcType=$bcArray[0];
	
	$bcCount++;
	print "$bcCount.|";
	print "$bcType|";
	
	print indexStat "$bcCount.|";
	print indexStat "$bcType|";
	
	#foreach $_ (@bcArray)
	#{
	#	print "$_|";
	#}
	#print "\n";
	
	if ($bcType !~ /^(barcode|dualBarcodes)$/i){ die "$bcType is not a valid barcode setting type!!\n"; }
	
	if ($bcType =~/^barcode/i){ $bNum=1; $bcArray[1]=uc $bcArray[1]; $bcCheck{$bcArray[1]}++;}
	elsif ($bcType =~/^dualBarcodes/i){ $bNum=2; $bcArray[1]=uc $bcArray[1]; $bcArray[2]=uc $bcArray[2]; $bcCheck{$bcArray[1].$bcArray[2]}++; }
	
	if ($itemNum < ($bNum+2)){die "$line\n$bNum/$itemNum Something missing! should be seqType & barcodes & sampleID!!\n\n";}
	elsif ($itemNum > ($bNum+2)){die "$line\n$bNum/$itemNum Got extra items ! should be seqType & barcodes & sampleID!!\n\n";}

	$barcodeTypeCount{$bcType}++;
	
	
	$sampleID=$bcArray[$bNum + 1]; if ($sampleID eq ""){$sampleID=$bcArray[$itemNum - 1];}
	if ($sampleID !~ /[A-Za-z0-9\-_]/){print "Warning: SampleID '$sampleID' should contain only alphanumeric characters & dash(-), not underline(_)! "; $sampleIDerr=1;}
	#else { push @sampleIDs,$sampleID; }
	
	
	for ($i=1;$i<= $bNum ;$i++)
	{
		#$bcArray[$i]=trim($bcArray[$i]);
		if ($bcArray[$i] =~ /[^ACGTUMRWSYKVHDBXN]/i){print "Warning: SampleID '$sampleID' barcode$i '$bcArray[$i]' contains invalid characters! "; $barcodeErr=1;}
		else
		{
			#$barcode{$sampleID}[$i]=$bcArray[$i];
			$barcodeLen[$i]=length($bcArray[$i]);
			$barcodeSeq{$i}{$bcArray[$i]}=$sampleID."_"."$i";
			
			print "$bcArray[$i]|$barcodeLen[$i]|";
			print indexStat "$bcArray[$i]|$barcodeLen[$i]|";
			#if ($i==1){$bcR1{$sampleID}=$bcArray[$i]; }#$bcR1Neg{$bcArray[0]}=cmplr(cmplc($bcArray[$i]));}  #R1 5'->3' at 5' => bcR1Pos; R1 5'->3' at 3' => bcR1Neg
			#if ($i==2){$bcR2{$sampleID}=$bcArray[$i]; }  #$bcR2Pos{$bcArray[0]}=$bcArray[$i];  #R2 5'->3' at 5' => bcR2Pos; R2 5'->3' at 3' => bcR2Neg;  need to convert bcR?Neg to reverse-complement to match reads
			
		
		}
	}
	
	if ($bcCheck{$bcArray[1].$bcArray[2]} >1)
	{
		print "  ==> Duplicated!!  ";
		print indexStat "  ==> Duplicated!!  ";
		$barcodeDupCheck=1;
	}

	
	if ($bNum==2){ $b2Sample{$bcArray[1]}{$bcArray[2]}=$sampleID; }
	else { $b2Sample{$bcArray[1]}=$sampleID; }
	
	print "$sampleID\n";
	print indexStat "$sampleID\n";
}
close (IN);
print indexStat "\n";

if ($sampleIDerr==1){die "There is at least one sampleID with invalid charactors!!\n\n";}
if ($barcodeErr==1){die "There is at least one barcode with invalid bases!!\n\n";}
if ($barcodeDupCheck==1)
{
	print "Duplicated barcode combination detected!! Are you sure to continue? (Y or N)"; $barcodeDupCheck=<STDIN>;
	if ($barcodeDupCheck !~ /[Yy]/){ die "Stop due to duplication of barcode cobinations!!\n\n"; }	
}


# put 5' barcode & 3' barcode into arrays
@Fbarcodes= keys %{$barcodeSeq{1}};
@Rbarcodes= keys %{$barcodeSeq{2}};


foreach $bType (sort keys %barcodeTypeCount)
{
	print "$bType: $barcodeTypeCount{$bType} found! with length: $barcodeLen[1] & $barcodeLen[2]\n";
	print indexStat "$bType: $barcodeTypeCount{$bType} found! with length: $barcodeLen[1] & $barcodeLen[2]\n";
}
print "\n";
print indexStat "\n";
	
for ($i=1;$i<= $Rnum ;$i++)
{
	$L[$i]=0;
	print "Reading $fq[$i]...\n";
	open (IN, "$fq[$i]") || die "$fq[$i] can not be openned.\n"; 
	while(defined($sT=<IN>))
	{ 
	  	#@NS500238:22:H03A2AFXX:1:11101:10176:1030 1:N:0:CGCTCATT
	  	#@M01199:14:000000000-AA4WM:1:1101:9448:1351 1:N:0:1
	  	#@HWI-ST615:245:D15FUACXX:5:1209:13303:100127 1:N:5122:GATCAG
	  	#@FCC011VACXX:2:1101:1385:1917#CGTAGGAC/1
			#@FCC011VACXX:2:1101:1385:1917#CGTAGGAC/2
#			if ($sT =~ /^@(.+)(\/[12]|\s+[12]:[YN]:([0-9]+):([0-9ATGCN]+))$/) 
			if ($sT =~ /^@(.+)(\/[12]|\s+[12]:[YN]:([0-9]+):(.+))$/)
			{
				$L[$i]++;
				$readLocation=$1; $readSuffix=$2;
				#$readLocation=~ s/:/-/g; $readLocation=~ s/_/-/g;  # replace : to -, and _ to -
				
				$readSuffix{$readLocation}[$i]=$readSuffix;  #print "$readLocation|$readSuffix{$readLocation}[$i]\n";

				if ($i==1){push @readLocationArray,$readLocation;}
				
				$sL = <IN>; chomp $sL;
				$readSeq{$readLocation}[$i]=$sL;
					
				$qT= <IN>; 
				if ($qT !~ /^\+/){print "invalid! @".$readLocation."'s quality title line is : $qT\nPress any key to continue...\n"; <STDIN>;}
				else
				{
					$qL= <IN>;  chomp $qL;
					$readQul{$readLocation}[$i]=$qL;
				}
		  }
		  else
		  {
		  	die "invalid! format not recognized!!\n$sT\n\n";
		  }
		  #@HWI-ST491:159:D0B5KACXX:7:1101:1209:2169 1:N:0:CGATGT
		  #@HWI-ST491:159:D0B5KACXX:7:1101:1209:2169 2:N:0:CGATGT
		  #@HWI-ST491:159:D0B5KACXX:2:1101:1070:2118 1:N:0:
	
		  
		  
		  #if ($L1==1){print "This fastq format is from CASAVA $fqType{$readLocation} !\n";}
	}
	close (IN);
	print "There are $L[$i] reads in $fq[$i]\n\n";
}


$indexOutputFile="indexOutput.$bFile";
open indexOut, ">$outputFoler/$indexOutputFile" or die $!;


$rCount=0; $udCount=0; $sampleRP=0; $undeterminedRP=0;
$orient{"F"}=$orient{"R"}=$orient{"NA"}=0;
foreach $readLocation (@readLocationArray)
{
	$rCount++;
	if ($bNum==2)
	{
		# if R1,R2 barcode length are the same, then ($R1headSeqF,$R2headSeqF) = ($R2headSeqR,$R1headSeqR)
		if ($barcodeLen[1]>0)
		{
			$R1headSeqF=uc substr($readSeq{$readLocation}[1],0,$barcodeLen[1]);
			$R2headSeqF=uc substr($readSeq{$readLocation}[2],0,$barcodeLen[1]);
		}
		else {$R1headSeqF="";$R2headSeqF="";}
		
		if ($barcodeLen[2]>0)
		{
			$R2headSeqR=uc substr($readSeq{$readLocation}[2],0,$barcodeLen[2]);
			$R1headSeqR=uc substr($readSeq{$readLocation}[1],0,$barcodeLen[2]);
		}
		else {$R2headSeqR="";$R1headSeqR="";}

		if (defined($b2Sample{$R1headSeqF}{$R2headSeqR}))
		{
			$sampleID=$b2Sample{$R1headSeqF}{$R2headSeqR}; $orient{"F"}++;
			
			$readSeqTrimBC{$readLocation}[1]=	substr($readSeq{$readLocation}[1],$barcodeLen[1]);
			$readSeqTrimBC{$readLocation}[2]=	substr($readSeq{$readLocation}[2],$barcodeLen[2]);
			
			$readQualTrimBC{$readLocation}[1]=	substr($readQul{$readLocation}[1],$barcodeLen[1]);
			$readQualTrimBC{$readLocation}[2]=	substr($readQul{$readLocation}[2],$barcodeLen[2]);
			
			print indexOut "$readLocation\t$sampleID\t$R1headSeqF\t$R2headSeqR\n";
		}
		elsif (defined($b2Sample{$R2headSeqF}{$R1headSeqR}))
		{
			$sampleID=$b2Sample{$R2headSeqF}{$R1headSeqR}; $orient{"R"}++;
			
			$readSeqTrimBC{$readLocation}[1]=	substr($readSeq{$readLocation}[2],$barcodeLen[1]);
			$readSeqTrimBC{$readLocation}[2]=	substr($readSeq{$readLocation}[1],$barcodeLen[2]);
			
			$readQualTrimBC{$readLocation}[1]=	substr($readQul{$readLocation}[2],$barcodeLen[1]);
			$readQualTrimBC{$readLocation}[2]=	substr($readQul{$readLocation}[1],$barcodeLen[2]);
			
			print indexOut "$readLocation\t$sampleID\t$R2headSeqF\t$R1headSeqR\n";
		}
		else 
		{ 
			#$R1Fmatch="";$R2Fmatch="";$R1Rmatch="";$R2Rmatch="";
			# ["initial_position=1","final_position=20","position_range=0"]
			# By setting the "position_range" to be zero you can limit (anchor) the operation to happen only once (if a match is possible) at the position.
			
			
			if ($barcodeLen[1]>0)
			{ 
				#@R1Fmatches = amatch($R1headSeqF, [ "S1","I1","D1" ], @Fbarcodes);
				
				if ($barcodeLen[1]>=16) { @R1Fmatches = amatch($R1headSeqF, [ "S2","I1","D1" ], @Fbarcodes); }
				elsif ($barcodeLen[1]>8) { @R1Fmatches = amatch($R1headSeqF, [ "S1","I1","D1" ], @Fbarcodes); }
				else { @R1Fmatches = amatch($R1headSeqF, [ "S1","I0","D0" ], @Fbarcodes); }
				
				$R1FmatchNum=scalar(@R1Fmatches); 
				if ($R1FmatchNum==0 || $R1FmatchNum >1){$R1Fmatches[0]="-NM-";} 
			}
			else {$R1Fmatches[0]="";}
			
			if ($barcodeLen[2]>0)
			{ 
				#@R2Rmatches = amatch($R2headSeqR, [ "S1","I1","D1" ], @Rbarcodes); 
				
				if ($barcodeLen[1]>=16) { @R2Rmatches = amatch($R2headSeqR, [ "S2","I1","D1" ], @Rbarcodes); }
				elsif ($barcodeLen[1]>8) { @R2Rmatches = amatch($R2headSeqR, [ "S1","I1","D1" ], @Rbarcodes); }
				else { @R2Rmatches = amatch($R2headSeqR, [ "S1","I0","D0" ], @Rbarcodes); }
				
				$R2RmatchNum=scalar(@R2Rmatches); 
				if ($R2RmatchNum==0 || $R2RmatchNum >1){$R2Rmatches[0]="-NM-";} 
			}
			else {$R2Rmatches[0]="";}
			
			#print "F:|$R1headSeqF|$R2headSeqF|=>|$R1Fmatches[0]|$R2Fmatches[0]|$R1FmatchNum|$R2FmatchNum\n";
			#if ($R1FmatchNum >1 || $R2FmatchNum >1){<STDIN>;}
			
			if ($barcodeLen[1]>0)
			{ 
				#@R2Fmatches = amatch($R2headSeqF, [ "S1","I1","D1" ], @Fbarcodes); 
				
				if ($barcodeLen[1]>=16) { @R2Fmatches = amatch($R2headSeqF, [ "S2","I1","D1" ], @Fbarcodes); }
				elsif ($barcodeLen[1]>8) { @R2Fmatches = amatch($R2headSeqF, [ "S1","I1","D1" ], @Fbarcodes); }
				else { @R2Fmatches = amatch($R2headSeqF, [ "S1","I0","D0" ], @Fbarcodes); }
					
				$R2FmatchNum=scalar(@R2Fmatches);
				if ($R2FmatchNum==0 || $R2FmatchNum >1){$R2Fmatches[0]="-NM-";} 
			}
			else {$R2Fmatches[0]="";}
			
			if ($barcodeLen[2]>0)
			{ 
				#@R1Rmatches = amatch($R1headSeqR, [ "S1","I1","D1" ], @Rbarcodes); 
				
				if ($barcodeLen[1]>=16) { @R1Rmatches = amatch($R1headSeqR, [ "S2","I1","D1" ], @Rbarcodes); }
				elsif ($barcodeLen[1]>8) { @R1Rmatches = amatch($R1headSeqR, [ "S1","I1","D1" ], @Rbarcodes); }
				else { @R1Rmatches = amatch($R1headSeqR, [ "S1","I0","D0" ], @Rbarcodes); }
				
				$R1RmatchNum=scalar(@R1Rmatches);
				if ($R1RmatchNum==0 || $R1RmatchNum >1){$R1Rmatches[0]="-NM-";} 
			}
			else {$R1Rmatches[0]="";}

			#print "R:|$R1headSeqR|$R2headSeqR|=>|$R1Rmatches[0]|$R2Rmatches[0]|$R1RmatchNum|$R2RmatchNum\n";
			#if ($R1RmatchNum >1 || $R2RmatchNum >1){<STDIN>;}
			
			if (defined($b2Sample{$R1Fmatches[0]}{$R2Rmatches[0]}))
			{
				$sampleID=$b2Sample{$R1Fmatches[0]}{$R2Rmatches[0]}; $orient{"F"}++;
				
				$readSeqTrimBC{$readLocation}[1]=	substr($readSeq{$readLocation}[1],$barcodeLen[1]);
				$readSeqTrimBC{$readLocation}[2]=	substr($readSeq{$readLocation}[2],$barcodeLen[2]);
				
				$readQualTrimBC{$readLocation}[1]=	substr($readQul{$readLocation}[1],$barcodeLen[1]);
				$readQualTrimBC{$readLocation}[2]=	substr($readQul{$readLocation}[2],$barcodeLen[2]);
				
				print indexOut "$readLocation\t$sampleID\t$R1Fmatches[0]\t$R2Rmatches[0]\t*\t$R1headSeqF\t$R2headSeqR\n";
			}
			elsif (defined($b2Sample{$R2Fmatches[0]}{$R1Rmatches[0]}))
			{
				$sampleID=$b2Sample{$R2Fmatches[0]}{$R1Rmatches[0]}; $orient{"R"}++;
			
				$readSeqTrimBC{$readLocation}[1]=	substr($readSeq{$readLocation}[2],$barcodeLen[1]);
				$readSeqTrimBC{$readLocation}[2]=	substr($readSeq{$readLocation}[1],$barcodeLen[2]);
				
				$readQualTrimBC{$readLocation}[1]=	substr($readQul{$readLocation}[2],$barcodeLen[1]);
				$readQualTrimBC{$readLocation}[2]=	substr($readQul{$readLocation}[1],$barcodeLen[2]);
				
				print indexOut "$readLocation\t$sampleID\t$R2Fmatches[0]\t$R1Rmatches[0]\t*\t$R2headSeqF\t$R1headSeqR\n";
			}
			else
			{
				$udCount++;
				$sampleID ="Undetermined"; 
				$orient{"NA"}++;
				
				if (defined($barcodeSeq{1}{$R1headSeqF})  &&  defined($barcodeSeq{1}{$R2headSeqF}) ){$sampleID.="FF"; }
				elsif (defined($barcodeSeq{2}{$R1headSeqR})  &&  defined($barcodeSeq{2}{$R2headSeqR}) ){$sampleID.="RR"; }
				elsif (  $R1Fmatches[0] ne "-NM-"  &&  $R2Rmatches[0] eq "-NM-" ){$sampleID.="1F"; $R1headSeqF=$barcodeSeq{1}{$R1Fmatches[0]}; }
				elsif (  $R1Fmatches[0] eq "-NM-"  &&  $R2Rmatches[0] ne "-NM-" ){$sampleID.="2R"; $R2headSeqR=$barcodeSeq{2}{$R2Rmatches[0]}; }
				elsif (  $R2Fmatches[0] ne "-NM-"  &&  $R1Rmatches[0] eq "-NM-" ){$sampleID.="2F"; $R2headSeqF=$barcodeSeq{1}{$R2Fmatches[0]}; }
				elsif (  $R2Fmatches[0] eq "-NM-"  &&  $R1Rmatches[0] ne "-NM-" ){$sampleID.="1R"; $R1headSeqR=$barcodeSeq{2}{$R1Rmatches[0]}; }
				
				#$R1headSeqT=uc substr($readSeq{$readLocation}[1],0,30);  
				#$R2headSeqT=uc substr($readSeq{$readLocation}[2],0,30);
				#print "$udCount/$rCount. $readLocation => $sampleID\t$R1headSeqF|$R2headSeqR\t$R2headSeqF|$R1headSeqR\t$R1headSeqT\t$R2headSeqT\n";
	
				$readSeqTrimBC{$readLocation}[1]=	$readSeq{$readLocation}[1];
				$readSeqTrimBC{$readLocation}[2]=	$readSeq{$readLocation}[2];
				
				$readQualTrimBC{$readLocation}[1]=	$readQul{$readLocation}[1];
				$readQualTrimBC{$readLocation}[2]=	$readQul{$readLocation}[2];
			}
		}
		
		if ($sampleID !~ /^Undetermined/){$sampleRP++;}
		else {$undeterminedRP++;}
		
		$countSampleIDnum{$sampleID}++;
		
		push @{$readsUnderSampleID{$sampleID}},$readLocation;
		
	}
	elsif ($bNum==1)
	{
		
		
		
		
		
	}
	else
	{
		die "$bNum  Something wrong! \n";
	}

}
close (indexOut);


#================
$sampleListFile=$bFile;
if ($sampleListFile =~ /barcodes/i){ $sampleListFile =~ s/barcodes/sampleFastqList/i ;} else {$sampleListFile="sampleFastqList.txt";}
open sampleListOut, ">$outputFoler/$sampleListFile" or die $!;


$SampleCount=0;
foreach $sampleID (sort keys %countSampleIDnum)
{
	$SampleCount++;
	$baseSpaceName1=$sampleID."_S".$SampleCount."_L001_R1_001.fastq";
	$baseSpaceName2=$sampleID."_S".$SampleCount."_L001_R2_001.fastq";
	
	$fq1OutFileName=$sampleID.".R1".$fqsuffix;
	$fq2OutFileName=$sampleID.".R2".$fqsuffix;
	
	print "Writing fastq files for $sampleID ...\n";

	open fq1Out, ">$outputFoler/$fq1OutFileName" or die $!;
	open fq2Out, ">$outputFoler/$fq2OutFileName" or die $!;
	
	foreach $readLocation (@{$readsUnderSampleID{$sampleID}})
	{
			print fq1Out "@".$readLocation.$readSuffix{$readLocation}[1]."\n";
			print fq1Out $readSeqTrimBC{$readLocation}[1]."\n";
			print fq1Out "+\n";
			print fq1Out $readQualTrimBC{$readLocation}[1]."\n";
		
			print fq2Out "@".$readLocation.$readSuffix{$readLocation}[2]."\n";
			print fq2Out $readSeqTrimBC{$readLocation}[2]."\n";
			print fq2Out "+\n";
			print fq2Out $readQualTrimBC{$readLocation}[2]."\n";
	}
	close fq1Out;
	close fq2Out;

	if (defined($outputBSFoler))
	{
		`cp -f $outputFoler/$fq1OutFileName $outputBSFoler/$baseSpaceName1`;
		`cp -f $outputFoler/$fq2OutFileName $outputBSFoler/$baseSpaceName2`;
		`gzip $outputBSFoler/$baseSpaceName1 &`;
		`gzip $outputBSFoler/$baseSpaceName2 &`;
	}

	print sampleListOut "$sampleID\t$fq1OutFileName\t$fq2OutFileName\n";
}
close (sampleListOut);




foreach $bType (sort keys %barcodeTypeCount)
{
	print "$bType: $barcodeTypeCount{$bType} found!\n";
	print indexStat "$bType: $barcodeTypeCount{$bType} found!\n";
	
}
print "\n";

print "sampleID\treadPairsFound\n";
print indexStat "sampleID\treadPairsFound\n";

$i=0;
foreach $sampleID (sort keys %countSampleIDnum)
{
	$i++;
	print "$i. $sampleID\t$countSampleIDnum{$sampleID}\n";
	print indexStat "$i\t$sampleID\t$countSampleIDnum{$sampleID}\n";
}
print "\n\n";
print indexStat "\n\n";

print "allSampleReadPairs\t$sampleRP\n";
print "allUndeterminedReadPairs\t$undeterminedRP\n\n";


print "Orientation\tCount:\n";
print "F\t".$orient{"F"}."\n";
print "R\t".$orient{"R"}."\n";
print "NA\t".$orient{"NA"}."\n";
print "\n\n";

print indexStat "allSampleReadPairs\t$sampleRP\n";
print indexStat "allUndeterminedReadPairs\t$undeterminedRP\n\n";


print indexStat "Orientation\tCount:\n";
print indexStat "F\t".$orient{"F"}."\n";
print indexStat "R\t".$orient{"R"}."\n";
print indexStat "NA\t".$orient{"NA"}."\n";
print indexStat "\n\n";



$date2=`date`;
print indexStat $date2;

close (indexStat);

print "\nDone!\n\n";

	if (defined($outputBSFoler))
	{
		print "Please check also folder '$outputBSFoler' for fastq files in standard Illumina format!!\n\n";
		
	}

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

