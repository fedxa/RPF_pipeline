#!/bin/sh
# tophat with -G options generates wrong strand XS:A:[+-] tags for antisence reads
# This causes problems for further invocation of cufflinks
# Fix them form the strand in the alignment flag
# https://genomebytes.wordpress.com/2013/07/08/fixing-the-xs-tag-in-tophat-output-bug-fixing/
#
# TODO: Should we add the tag to the sequences lacking one? (Seems not important?)
#       Add statistics dump to the stderr

## We use NTHREADS-1 compression threads, by default 3 total CPU assumed
NTHREADS=${NTHREADS:=3}
NTHREADS1=$[ $NTHREADS - 1 ]

samtools view -h $1 \
    | awk '{ if($2%256==16) {sub(/XS\:A\:\+/,"XS:A:-",$0); print $0};
             if($2%256==0)  {sub(/XS\:A\:\-/,"XS:A:+",$0); print $0}; }' \
    | samtools view -@$NTHREADS1 -h -S -b -
