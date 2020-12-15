#!/usr/bin/perl -w
#
# Process VEP output (VCF)
# /export/EC1680U/Daniel_WES/GNND-Zheng160127/DNA15-064/output/vqsr
# variant_effect_predictor.pl -i DNA15-064.filtered.vcf --cache --everything --dir /export/arrayPRO1/commonSoftwares/ensembl-tools-release-81/.vep/ --cache_version 81 --pick --vcf --fork 30 --ASSEMBLY GRCh38
# ./vcf2txt.v2.pl variant_effect_output.txt > DNA15-064.variant_effect_output.txt

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
      #print "$tmp\n\n";	  
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
      $finalfield = "#CHROM\tPOS\tID\tREF\tALT\tSET\tQUAL\tFILTER\tZygosity\tAlleleFrequency\tDepth\t".$field; 
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
    $info=~/.*(set=.*);.*/;
    $set = $1;    
    
    @semicolon = split(/;/, $info); 
    for ($i=0; $i< scalar @semicolon -1; $i++) {
      if ($semicolon[$i]=~/AF=(.*)/) {
        $af = $1;          # allele frequency        
      }
      if ($semicolon[$i]=~/DP=(.*)/) {
        $depth =  $1;         # depth	    
      }	  
    }
    $tmp   = $semicolon[-1];
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
	
	$zygo = ""; 
	for ($idx=9; $idx < scalar @tab; $idx++) {               # sample
	  $samp = $tab[$idx];      
      @colon = split(/:/, $samp);
	  if ($colon[0] eq "./.") { next;}
      if ($colon[0] eq "1/1") {
        $zygo = $zygo."Hom|";
      } elsif ($colon[0] eq "0/1") {
        $zygo = $zygo."Het|";
      } else {
	    $zygo = $zygo.".|";
	  }
	}
	chop($zygo);
	
    print "$chr\t$pos\t$id\t$ref\t$alt\t$set\t$qual\t$filt\t$zygo\t$af\t$depth\t";
    for ($idx2=1; $idx2 < scalar @fields; $idx2++) {
      $output = $data->{$fields[$idx2]};		
      if ($fields[$idx2] eq "Condel") {		    	
        $newoutput = processvalue($output);				
      } elsif ($fields[$idx2] eq "SIFT") {	
	    #print "SIFT ".$output."AA\n";
	    $newoutput = processvalue($output);		
      } elsif ($fields[$idx2] eq "PolyPhen") {	
	    #print "PolyPhen ".$output."AA\n";
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
