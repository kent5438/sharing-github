#!/usr/bin/perl
use strict;
use warnings;
use File::Path qw(make_path);

my $input_taxa	= defined $ARGV[0]? $ARGV[0] : '';
my $mapping_file= defined $ARGV[1]? $ARGV[1] : '';
my $group_title	= defined $ARGV[2]? $ARGV[2] : 'Description';
my $output_folder=defined $ARGV[3]? $ARGV[3] : 'LEfSe';

die "perl $0 [input taxa txt] [mapping txt] ([group title: Description]) ([output folder: LEfSe])\n" if(!-f $input_taxa || !-f $mapping_file);
make_path($output_folder);


# read mapping file for group
my %group;
my %group_count;
my @sample;
{
my $idIndex = 0;
my $groupIndex = 0;
my %sampleCheck;
open(IN,'<'.$mapping_file);
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  my @array = split("\t",$_);

  if (@array < 4) {
    next;
  }

  if ($_ =~ /^#/) {
    for (my $i=0;$i<@array;$i++) {
      if ($array[$i] eq '#SampleID') {
        $idIndex = $i;
      }
      elsif ($array[$i] eq $group_title) {
        $groupIndex = $i;
      }
    }
  }
  else {
    $array[$groupIndex] =~ s/\s+/_/g; #replace space

    if ($sampleCheck{$array[$idIndex]} == 1) { next; }

    push(@sample,$array[$idIndex]);

    $group{$array[$idIndex]} = $array[$groupIndex];
    $group_count{$array[$groupIndex]} ++;

    $sampleCheck{$array[$idIndex]} = 1;
  }
}
close(IN);
}

die "LEfSe need two groups at least.\n" if(keys %group_count < 2);

# read taxa file
my %taxa;
{
open(IN,'<'.$input_taxa);
my @header;
while (<IN>) {
  chomp;
  if ($_ eq '') {
    next;
  }

  if ($_ =~ /^Taxon/) {
    @header = split("\t",$_);
    next;
  }

  my @array = split("\t",$_);

  if ((@header == 0) || ($#header != $#array)) {
  }

  $array[0] =~ s/\s//g;
  $array[0] =~ s/;/\|/g;

  if ($array[0] =~ /[kpcofgs]__$/) { next; }

  $array[0] =~ s/^[kpcofgs]__//g;
  $array[0] =~ s/\|[kpcofgs]__/\|/g;

  for (my $i=1;$i<@array;$i++) {
    if (!defined($group{$header[$i]}) || ($group{$header[$i]} eq '')) {
      next;
    }

    $taxa{$array[0]}{$header[$i]} += $array[$i];
  }
}
close(IN);
}

# output
open(OUT,'>'.$output_folder.'/LEfSe_input.txt');
print OUT join("\t",'SampleID', @sample)."\n";
print OUT join("\t",$group_title, @{\%group}{@sample})."\n";
foreach my $tax (sort{$a cmp $b} (keys %taxa)) {
  print OUT $tax;

  foreach my $sampleName (@sample) {
    my $val = $taxa{$tax}{$sampleName} + 0;

    print OUT "\t".$val;
  }

  print OUT "\n";
}
close(OUT);

# run LEfSe
`format_input.py $output_folder/LEfSe_input.txt $output_folder/LEfSe_formated.txt -u 1 -c 2 -o 1000000`;
my $log = `run_lefse.py $output_folder/LEfSe_formated.txt $output_folder/LEfSe_output.txt`;
die "No features with significant differences between the two classes\n" if($log =~ /No features with significant differences between the two classes/);
`plot_res.py $output_folder/LEfSe_output.txt $output_folder/LEfSe.png --width 15 --feature_font_size 10`;
`plot_cladogram.py $output_folder/LEfSe_output.txt  $output_folder/LEfSe_cladogram.png --format png --dpi 150  --left_space_prop 0`;
