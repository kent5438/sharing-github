#!/usr/bin/perl
use strict;
use warnings;

our $mapping_file  = defined $ARGV[0]? shift @ARGV : '';
our $phred_offset  = defined $ARGV[0]? shift @ARGV : 33; # 33 or 64
our $quality       = defined $ARGV[0]? shift @ARGV : 19; # mean Q20

# FLASH options
our $min_overlap   = defined $ARGV[0]? shift @ARGV : 10;        # -m,  The minimum required overlap length between two reads to provide a confident overlap.  Default: 10bp.
our $max_overlap   = defined $ARGV[0]? shift @ARGV : 250;        # -M,  Maximum overlap length expected in approximately 90% of read pairs.  It is by default set to 65bp, which works well for 100bp reads generated from a 180bp library.
our $mismatch_ratio= defined $ARGV[0]? shift @ARGV : 0.25;      # -x,  Maximum allowed ratio between the number of mismatched base pairs and the overlap length. Default: 0.25.

# AlienTrimmer options
our $trim_fasta    = defined $ARGV[0]? shift @ARGV : '';        # -b,  for single-end & pair-end
our $trim_fasta_f  = defined $ARGV[0]? shift @ARGV : '';        # -g, for pair-end
our $trim_fasta_r  = defined $ARGV[0]? shift @ARGV : '';        # -a, for pair-end
our $kmer          = defined $ARGV[0]? shift @ARGV : 10;        # -k,  k-mer decomposition, between 5 and 15
our $mismatches    = defined $ARGV[0]? shift @ARGV : ($kmer/2); # -m,  maximum mismatches, between 0 and 15
our $read_length   = defined $ARGV[0]? shift @ARGV : 150;        # -l,  minimum read length
our $phred_quality = defined $ARGV[0]? shift @ARGV : 20;        # -q,  Phred quality score cut-off, between 0 and 40
our $percentage    = defined $ARGV[0]? shift @ARGV : 0;         # -p,  minimum allowed percentage of correctly called nucleotides, between 0 and 100

our $thread        = defined $ARGV[0]? shift @ARGV : 32;
our $printGap      = defined $ARGV[0]? shift @ARGV : '';

# check parameters
&Useage();

# validate mapping file
my @cache = split(/\n/, &ValidateMappingFile());

# setup trimming tools
my $jar = `echo \$AlienTrimmer`; chomp($jar);
$jar = (-e $jar)? $jar : '/app_tools/AlienTrimmer.jar'; # set default path

# parse info in mapping file
my $idIndex = 0;
my $groupIndex = 0;
my %group = ();
my @fastqFiles = ();
foreach(@cache){
	my @column = split(/\t/, $_);
	if($_ =~ /^#(.+)/){
		@column = split(/\t/, $1);
		for my $i(0..$#column){
			$idIndex = $i if($column[$i] eq 'SampleID');
			$groupIndex = $i if($column[$i] eq 'FilePath');
		}
	}else{
		$group{$column[$idIndex]} = $column[$groupIndex];
		push @fastqFiles, split(/,/, $column[$groupIndex]);
	}
}

# B.1 QC
&PrintProcess('B.1 QC');
`rm -rf B.1_QC && mkdir B.1_QC`;

&PrintProcess('    fastaQC-rawdata');
my $cmd = "fastqc -o B.1_QC/fastqc-rawdata -f fastq -q -t $thread ". join(" ", @fastqFiles);
`mkdir B.1_QC/fastqc-rawdata && $cmd`;
#&PrintProcess('    multiQC-rawdata');
#`multiqc -f -q -o B.1_QC/multiqc-rawdata B.1_QC/fastqc-rawdata`;

# trimming
&PrintProcess('    trimming');
if(-e $jar && (-e $trim_fasta || -e $trim_fasta_f)){
    my $trimmed_PE = 0;
	my $trimmed_SE = 0;
	`mkdir B.1_QC/trimming`;
	my $pwd = `pwd`;
	chomp($pwd);
	
	my @trimmedFastq = ();
	foreach my $sample (keys %group){
		my @fastq = split(/,/, $group{$sample});
		my $generalOption = "cutadapt ";
        
		# (MiSeq) Pair-end
		if($#fastq > 0){
            my $c_option = -e $trim_fasta ? "-m $read_length --discard-untrimmed -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC" : "";
            if(-e $trim_fasta_f){
                $c_option .= -e $trim_fasta_f ? " -g $trim_fasta_f" : "";
                $c_option .= -e $trim_fasta_r ? " -a $trim_fasta_r" : "";
            }
            
			$fastq[0] =~ /\/([^\/]+)$/; my $f_name = $1;
			$fastq[1] =~ /\/([^\/]+)$/; my $r_name = $1;
			`$generalOption $c_option -o B.1_QC/trimming/$f_name -p B.1_QC/trimming/$r_name $fastq[0] $fastq[1] > /dev/null`;
            $group{$sample} = "$pwd/B.1_QC/trimming/$f_name,$pwd/B.1_QC/trimming/$r_name";
			push @trimmedFastq, "$pwd/B.1_QC/trimming/$f_name", "$pwd/B.1_QC/trimming/$r_name";
			$trimmed_PE = 1;
		# (Pacbio) Single-end
		}else{
            last if !-e $trim_fasta; # single-end must have $trim_fasta
			$fastq[0] =~ /\/([^\/]+)$/; my $f_name = $1;
			`$generalOption -m 1000 -g GCAGTCGAACATGTAGCTGACTCAGGTCACAGAGTTGATCMTGGCTCAG -g TGGATCACTTGTGCAAGCATCACATCGTAGGGTTACCTTGTTACGACTT -o B.1_QC/trimming/$f_name $fastq[0] > /dev/null`;
			$group{$sample} = "$pwd/B.1_QC/trimming/$f_name";
			push @trimmedFastq, "$pwd/B.1_QC/trimming/$f_name";
			$trimmed_SE = 1;
		}
#        $trimmed = 1;
	}
	
    if($trimmed_PE){
		&PrintProcess('    multiQC-rawdata');
		`multiqc -f -q -c /export/EC1680U/perl/bin/16S/MiSeq-16S_multiqc_config.yaml -o B.1_QC/multiqc-rawdata B.1_QC/fastqc-rawdata`;
        &PrintProcess('    fastaQC-trimmed');
        my $cmd = "fastqc -o B.1_QC/fastqc-trimmed -f fastq -q -t $thread ". join(' ', @trimmedFastq);
        `mkdir B.1_QC/fastqc-trimmed && $cmd`;
        &PrintProcess('    multiQC-trimmed');
        `multiqc -f -q -c /export/EC1680U/perl/bin/16S/MiSeq-16S_multiqc_config.yaml -o B.1_QC/multiqc-trimmed B.1_QC/fastqc-trimmed`;
    } elsif($trimmed_SE){
		&PrintProcess('    multiQC-rawdata');
		`multiqc -f -q -c /export/EC1680U/perl/bin/16S/PB-16S_multiqc_config.yaml -o B.1_QC/multiqc-rawdata B.1_QC/fastqc-rawdata`;
		&PrintProcess('    fastaQC-trimmed');
		my $cmd = "fastqc -o B.1_QC/fastqc-trimmed -f fastq -q -t $thread ". join(' ', @trimmedFastq);
		`mkdir B.1_QC/fastqc-trimmed && $cmd`;
		&PrintProcess('    multiQC-trimmed');
		`multiqc -f -q -c /export/EC1680U/perl/bin/16S/PB-16S_multiqc_config.yaml -o B.1_QC/multiqc-trimmed B.1_QC/fastqc-trimmed`;
	}
}

# B.2 join paried-end reads
&PrintProcess('B.2 join paried-end reads');
`rm -rf B.2_join_fastq && mkdir B.2_join_fastq`;
my %case = (); # for B.3
foreach my $sample (keys %group){
	my @fastq = split(/,/, $group{$sample});
	if($#fastq > 0){
		`flash $fastq[0] $fastq[1] -p $phred_offset -m $min_overlap -M $max_overlap -x $mismatch_ratio -d B.2_join_fastq/$sample -o $sample -t $thread > B.2_join_fastq/$sample.log 2>&1`;
		$case{"B.2_join_fastq/$sample/$sample.extendedFrags.fastq"} = $sample;
	}else{
		$case{$group{$sample}} = $sample;
	}
}

# B.3 demultiplex and quality filtering
&PrintProcess('B.3 demultiplex and quality filtering');
`rm -rf B.3_demultiplex_quilty_filter`;
my $inputFiles = '';
my $sampleIds  = '';
foreach my $fileName (keys %case){
	$inputFiles .= "$fileName,";
	$sampleIds  .= "$case{$fileName},";
}
$inputFiles =~ s/,$//g; $sampleIds =~ s/,$//g;
#`split_libraries_fastq.py --phred_quality_threshold $quality --phred_offset $phred_offset --barcode_type 'not-barcoded' -i $inputFiles --sample_ids $sampleIds -o B.3_demultiplex_quilty_filter`;
system('perl /export/EC1680U/perl/bin/16S/mergeSeq.pl '.$mapping_file.' B.3_demultiplex_quilty_filter');

sub Useage(){
    my $checker = 0;
    $checker = $checker || !-f $mapping_file;
    $checker = $checker || ($phred_offset != 33 && $phred_offset != 64);
    $checker = $checker || $quality !~ /^[0-9]+$/;
    
    $checker = $checker || $min_overlap !~ /^[0-9]+$/;
    $checker = $checker || $max_overlap !~ /^[0-9]+$/;
    $checker = $checker || $mismatch_ratio !~ /^[0-9.]+$/ || $mismatch_ratio > 1; # -x <= 1
    
    $checker = $checker || $kmer !~ /^[0-9]+$/ || $kmer < 5 || $kmer > 15; # -k[5-15]
    $checker = $checker || $mismatches !~ /^[0-9]+$/ || $mismatches < 0 || $mismatches > 15; # -m[0-15]
    $checker = $checker || $read_length !~ /^[0-9]+$/;
    $checker = $checker || $phred_quality !~ /^[0-9]+$/ || $phred_quality < 0 || $phred_quality > 40; # -q[0-40]
    $checker = $checker || $percentage !~ /^[0-9]+$/ || $percentage < 0 || $percentage > 100; # -p[0-100]
    
    $checker = $checker || $thread !~ /^[0-9]+$/;
    
    my $usage = "perl $0";
    $usage .= " [mapping file]";
    $usage .= " ([phred_offset: 33])";
    $usage .= " ([quality: 19])";    
    $usage .= " ([min overlap: 10])";
    $usage .= " ([max overlap: 65])";
    $usage .= " ([mismatch ratio: 0.25])";    
    $usage .= " [trim fasta: '']";
    $usage .= " [forware trim fasta: '']";
    $usage .= " [reverse trim fasta: '']";
    $usage .= " ([kmer: 10])";
    $usage .= " ([mismatches: kmer/2])";
    $usage .= " ([read length: 15])";
    $usage .= " ([phred quality: 20])";
    $usage .= " ([allowed percentage: 0])";
    $usage .= " ([thread: 8])";
    $usage .= " ([log gap: ''])";
    
    die "$usage\n" if($checker);
}

sub ValidateMappingFile(){
	my $cache   = '';
	
	my $infoCheck = 0;
	my $idIndex = 0;
	my $fileIndex = 0;
	my %sampleID = ();
	my $totalColumn = 0;
	open MappingFile, $mapping_file;
		while(<MappingFile>){
			chomp;
			my @column = split(/\t/, $_, -1);
			next if($#column < 4); # at least: SampleID, BarcodeSequence, LinkerPrimerSequence, Description, FilePath
			$cache .= "$_\n";
			if($_ =~ /^#(.+)/){
				@column = split(/\t/, $1);
				$totalColumn = $#column;
				my $titleCheck = 0;
				for my $i(0..$#column){				
					if($column[$i] eq 'SampleID'){
						$idIndex = $i;
						$titleCheck += 1;
					}elsif($column[$i] eq 'BarcodeSequence'){
						$titleCheck += 10;
					}elsif($column[$i] eq 'LinkerPrimerSequence'){
						$titleCheck += 100;
					}elsif($column[$i] eq 'Description'){
						$titleCheck += 1000;
					}elsif($column[$i] eq 'FilePath'){
						$fileIndex = $i;
						$titleCheck += 10000;
					}
				}
				die "MappingFile must include columns: 'SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'Description', 'FilePath'.\n" if($titleCheck != 11111);
				$infoCheck = 1;
			}else{
				die "Duplicate SampleID.\n" if(defined($sampleID{$column[$idIndex]}));
				$sampleID{$column[$idIndex]} = 1;
				
				die "Inconsistent column number occured." if($#column != $totalColumn);
				
				die "SampleID column must not contain '.' or '|'." if($column[$idIndex] =~ /[.|]/);
				$infoCheck = 2 if($infoCheck);
				
				my @files = split(/,/, $sampleID{$column[$idIndex]});
				die "Too much file, pair-end include two files, single-end include one file.\n" if($#files >1);
				die "FilePath column must not be empty.\n" if($#files < 0);
			}
		}
	close MappingFile;	
	die "Empty mapping file\n" if($infoCheck != 2);
	
	return $cache;
}

sub PrintProcess(){
	my $str = shift @_;
	print "$printGap$str\n";
}
