#!/bin/sh
# Simple an dirty multihit randomization
# Select single hits (based on NH:i:# tag)
# Select one of the double hits (at 0.5 probability)
# Drop reads with >2 alignments

if [ $# -lt 2 ]; then
    echo "Usage: $0 in.bam out.bam"
    exit 1
fi
in=$1
out=$2

basedir=`dirname $0`

(
    samtools view -H $in
    samtools view $in |grep -E 'NH:i:1([^0-9]|$)'
    samtools view $in |grep -E 'NH:i:2([^0-9]|$)' | sort | python3 $basedir/randomizeDoubleHits.py
) | samtools view -b > $in.tmp.bam
samtools sort $in.tmp.bam > $out
rm $in.tmp.bam
