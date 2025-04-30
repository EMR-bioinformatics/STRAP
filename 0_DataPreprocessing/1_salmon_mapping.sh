#!/bin/bash

outdir=$1
threads=$2
fastaPath=$3
sampleList=$4

mkdir $outdir

while read filename; do

/tmp/salmon-latest_linux_x86_64/bin/salmon quant -p $threads -i "/tmp/gencodetranscriptsindexv29mod/" \
 --gcBias --libType A -1 $fastaPath/$filename"_R1_001.fastq.gz" -2 $fastaPath/$filename"_R2_001.fastq.gz" -o  $outdir/$filename/"salmon"

done < $sampleList

