#!/bin/sh
# tophat with -G options generates wrong strand XS:A:[+-] tags
# Fix them form the strand in the alignment flag
# https://genomebytes.wordpress.com/2013/07/08/fixing-the-xs-tag-in-tophat-output-bug-fixing/
#
# TODO: Should we add the tag to the sequences lacking one?
#       Seems not important?

## We use NTHREADS-1 compression threads, by default 3 total CPU assumed
NTHREADS=${NTHREADS:=3}
NTHREADS1=$[ $NTHREADS - 1 ]

samtools view -h $1 \
    | awk '{ if($2%256==16) {sub(/XS\:A\:\+/,"XS:A:-",$0); print $0};
             if($2%256==0)  {sub(/XS\:A\:\-/,"XS:A:+",$0); print $0}; }' \
    | samtools view -@$NTHREADS1 -h -S -b -


# samtools view -h $1 \
#     | awk '/XS:A:+/ { if($2%256==16) {sub(/XS\:A\:\+/,"XS:A:-",$0); wrongap++} else {goodap++}; print; next }
#            /XS:A:-/ { if($2%256==0)  {sub(/XS\:A\:\-/,"XS:A:+",$0); wrongam++} else {goodam++}; print; next }
#            /^#/        { print; next }
#                        { print; noa++ }
#            END {print "wrongap:",wrongap,"\nwrongam:",wrongam,"\ngoodap:",goodap,"\ngoodam:",goodam,"\nnoa:",noa > "/dev/stderr"}' \
#     | samtools view -@$NTHREADS1 -h -S -b -
