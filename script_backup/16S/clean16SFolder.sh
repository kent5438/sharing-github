folder=$1
if [ ! $folder ]; then
  echo No argument;
  exit
fi
if [ ! -e $folder ]; then 
  echo $folder not exist;
  exit;
fi

cd $folder
#cd OTU_process
/bin/rm -fR 1_joinPairs 2_removePrimers 3_screenQuality 4_checkChimera
mkdir 1_joinPairs 2_removePrimers 3_screenQuality 4_checkChimera
ls 5_assignTaxa/*/* | grep -v 'biom$' | grep -v 'fasta$' | grep -v 'names$' | grep -v 'taxonomy$' | perl -e 'print "/bin/rm -f";while (<>) {chomp; print " ".$_;}print "\n";' | sh -
cd ../
du -hs *.16S | sort -h
cd ../
