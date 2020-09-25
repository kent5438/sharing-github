#! /bin/bash

asm=$1
if [ -z $asm ]; then
    echo "### ERROR: 10x assembly output is required!"
    exit 1
fi

supernova mkoutput \
    --asmdir="${asm}/outs/assembly" \
    --outprefix="${asm}/raw" \
    --style=raw

## megabubbles
supernova mkoutput \
    --asmdir="${asm}/outs/assembly" \
    --outprefix="${asm}/megabubbles" \
    --style=megabubbles

## pseudohap
supernova mkoutput \
    --asmdir="${asm}/outs/assembly" \
    --outprefix="${asm}/pseudohap" \
    --style=pseudohap

## pseudohap2
supernova mkoutput \
    --asmdir="${asm}/outs/assembly" \
    --outprefix="${asm}/pseudohap2" \
    --style=pseudohap2
