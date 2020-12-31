#! /bin/bash
### Run on Snake only (SSD + large local space)


echo "10x supernova has been under maintain!"

prefix=$1
if [ -z $prefix ]; then
	echo "### ERROR: 10x project ID is required!"
	exit 1
fi

# mkfastq (need to reformat samplesheet in advance)
#supernova mkfastq \
#	--run=/mnt/NFS/HiSeq4000/171128_K00236_0177_BHLGCNBBXX/ \
#	--localcores=32 \
#	--localmem=150G \
#	--lanes=8 \
#	--output-dir=mkfastq \
#	--barcode-mismatches=1 \
#	--csv /mnt/NFS/HiSeq4000/171128_K00236_0177_BHLGCNBBXX/tiny-bcl-simple-2.1.0.csv \
#	--ignore-dual-index

# run
supernova run \
	--id "${prefix}_supernova" \
	--fastqs="${prefix}/" \
	--localcores=32 \
	--localmem=250 \
	--maxreads=533000000 # The maxreads depend on the 56x/70x coverage and the read length


# mkoutput
## raw
supernova mkoutput \
	--asmdir="${prefix}_supernova/outs/assembly" \
	--outprefix="${prefix}_supernova/raw" \
	--style=raw

## megabubbles
supernova mkoutput \
	--asmdir="${prefix}_supernova/outs/assembly" \
	--outprefix="${prefix}_supernova/megabubbles" \
	--style=megabubbles

## pseudohap
supernova mkoutput \
	--asmdir="${prefix}_supernova/outs/assembly" \
	--outprefix="${prefix}_supernova/pseudohap" \
	--style=pseudohap

## pseudohap2
supernova mkoutput \
	--asmdir="${prefix}_supernova/outs/assembly" \
	--outprefix="${prefix}_supernova/pseudohap2" \
	--style=pseudohap2

