#!/usr/bin/perl -w
#
# Process VEP output (VCF)
# variant_effect_predictor.pl -i output/mutect/my.accept.vcf --cache --everything --dir /export/arrayPRO1/commonSoftwares/ensembl-tools-release-81/.vep/ --cache_version 81 --pick --vcf --fork 30 --offline --ASSEMBLY GRCh37 -o output/vep/variant_effect_output.txt
# ./vcf2txt.mutect.v1.pl variant_effect_output.txt > DNA15-064.variant_effect_output.txt

my $fn = shift;
my $header = "";
my $field = "";
my $finalfield = "";
my @fields;
my $idx  = 1;
my $idx2 = 1;
my $idx3 = 1;
my $check = 1;
my $tmp_value = "";
my $newoutput;

open (IN, "<$fn") or die "Can't open $fn file";
while (<IN>) {
  $ln = $_;
  chomp($ln);
  if ($ln=~/#/) {
    if ($ln=~/#CHROM/) {
      next;
    }
    $header = $header.$ln."\n";
    if ($ln=~/##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: (.*)\">/) {
      $tmp = $1;
      @fields = split(/\|/, $tmp);
      for ($idx=1; $idx < scalar @fields; $idx++) {
        if ($fields[$idx] eq "Condel") {
	      $field =$field."Condel_Desc\tCondel_value\t";
	    } elsif ($fields[$idx] eq "SIFT") {
	      $field =$field."SIFT_Desc\tSIFT_value\t";
	    } elsif ($fields[$idx] eq "PolyPhen") {
	      $field =$field."PolyPhen_Desc\tPolyPhen_value\t";
	    } elsif ($fields[$idx] eq "GMAF") {
    	  $field =$field."GMAF_Allele\tGMAF_Freq\t";
	    } elsif ($fields[$idx] eq "AFR_MAF") {
	      $field =$field."AFR_MAF_Allele\tAFR_MAF_Freq\t";
	    } elsif ($fields[$idx] eq "AMR_MAF") {
	      $field =$field."AMR_MAF_Allele\tAMR_MAF_Freq\t";
	    } elsif ($fields[$idx] eq "EAS_MAF") {
	      $field =$field."EAS_MAF_Allele\tEAS_MAF_Freq\t";
        } elsif ($fields[$idx] eq "EUR_MAF") {
	      $field =$field."EUR_MAF_Allele\tEUR_MAF_Freq\t";		  
	    } elsif ($fields[$idx] eq "SAS_MAF") {
	      $field =$field."SAS_MAF_Allele\tSAS_MAF_Freq\t";		  
	    } elsif ($fields[$idx] eq "AA_MAF") {
	      $field =$field."AA_MAF_Allele\tAA_MAF_Freq\t";		  
	    } elsif ($fields[$idx] eq "EA_MAF") {
	      $field =$field."EA_MAF_Allele\tEA_MAF_Freq\t";		  
	    } else {
          $field =$field.$fields[$idx]."\t";
	    }
      }
      chop($field);
      $field = $field."\n";	  
    }
  } else {
    if ($check ==1) {	
      $finalfield = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tZygosity\tAlleleFrequency\tDepth\t".$field; 
      print "$header$finalfield";
      $check = $check + 1;
    }
    @tab  = split(/\t/, $ln); 
    $chr  = $tab[0];                 # chromosome
    $pos  = $tab[1];                 # position
    $id   = $tab[2];                 # ID
    $ref  = $tab[3];                 # Reference Allele
    $alt  = $tab[4];                 # Alternative Allele
    $qual = $tab[5];                 # Quality
    $filt = $tab[6];                 # filter
    $info = $tab[7];                 # Info
    #$form = $tab[8];                # Format    
    $samp = $tab[9];                 # sample [Tumor]
        
    @colon = split(/:/, $samp); 
	$gt = $colon[0];                 # Genotype
	$ad = $colon[1];                 # Allelic depths for the ref and alt alleles in the order listed	
	@tmp = split(/,/, $ad);
	$depth = $tmp[0] + $tmp[1];      # depth
	$af = $colon[2];                 # Allele fraction of the event in the tumor
	
    
	@semicolon = split(/;/, $info);
	#for ($idx=0; $idx < scalar @semicolon; $idx ++) {
	#  if ($semicolo[$idx]=~/CSQ/) {
	$tmp = $semicolon[-1];    
	#  }
	#}
    
    $tmp =~/CSQ=(.*)/;
    $csq = $1;
    $csq =~s/\|/ \|/g;               # Very Important
    @pipe = split(/\|/, $csq);	  
    for ($idx2=1; $idx2 < scalar @fields; $idx2++) {	    
      $tmp_value = $pipe[$idx2];	        
      unless (defined $pipe[$idx2]) {
        $tmp_value = " ";
      }
      $data->{$fields[$idx2]} = $tmp_value;	    		
      #print $fields[$idx2]."\t".$tmp_value."\n"; # Important check point
    }	  
    
    if ($gt eq "1/1") {
      $zygo = "Hom";
    } 
    if ($gt eq "0/1") {
      $zygo = "Het";
    }
	
    print "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filt\t$zygo\t$af\t$depth\t";
    for ($idx2=1; $idx2 < scalar @fields; $idx2++) {
      $output = $data->{$fields[$idx2]};		
      if ($fields[$idx2] eq "Condel") {		    	
        $newoutput = processvalue($output);				
      } elsif ($fields[$idx2] eq "SIFT") {	
	$newoutput = processvalue($output);		
      } elsif ($fields[$idx2] eq "PolyPhen") {	
	$newoutput = processvalue($output);	
      } elsif ($fields[$idx2] eq "GMAF") {	
	$newoutput = processvalue($output);			
      } elsif ($fields[$idx2] eq "AFR_MAF") {	
	$newoutput = processvalue($output);							
      } elsif ($fields[$idx2] eq "AMR_MAF") {	
	$newoutput = processvalue($output);					
      } elsif ($fields[$idx2] eq "EAS_MAF") {	
	$newoutput = processvalue($output);						
      } elsif ($fields[$idx2] eq "EUR_MAF") {	
	$newoutput = processvalue($output);						
      } elsif ($fields[$idx2] eq "SAS_MAF") {	
	$newoutput = processvalue($output);							
      } elsif ($fields[$idx2] eq "AA_MAF") {	
	$newoutput = processvalue($output);						
      } elsif ($fields[$idx2] eq "EA_MAF") {	
	$newoutput = processvalue($output);							
      } else {	    
        $newoutput = $output."\t";	
      }	  
      print "$newoutput";
    }
    print "\n";
  } 
}
close IN;
sub processvalue {
  my $value = shift;
  my $newvalue;
    
  #if ($value=~/(\D+)\((.*)\)/) {
  if ($value=~/(\D+)\((.*)\)/) {
    $newvalue = "$1\t$2\t";
  } elsif ($value=~/(\D)\:(.*)/) {
    $newvalue = "$1\t$2\t";
  } else {
    $newvalue = "\t\t";
  }

  return $newvalue;
}

1;
