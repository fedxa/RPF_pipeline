#!/bin/bash
# Starts cutadapt on cluster in parallel

module add UHTS/Quality_control/cutadapt/1.8

if [ $# -lt 4 ]; then
    echo "Usage: $0 in_1.fastq in_2.fastq out_1.fatsq out_2.fastq [cutadapt-options]"
    echo ' Runs cutadapter on LSB cluster guessing properly number of cores/nodes'
    echo ' All files can be fastq, fastq.gz, fastq.bz2, fastq.xz (all should have the same type'
    echo ' cutadapt-options are various options (i.e. -a ADAPTER) transered to cutadapt without change'
    echo ' (use cutadapt --help to see the list)'
    echo ' Number of processes is deduced automatically under bsub, or taken from environment NCHUNKS'
    echo ' Use empty files "" for the second pair for single end mode'
    exit 1
fi
in1=$1
shift
in2=$1
shift
out1=$1
shift
out2=$1
shift

# Magic to guess number of hosts
if [ "x$LSB_HOSTS" = "x" ]; then
    NO_LSB=True
    if [ "x$NCHUNKS"  = "x" ]; then
	NCHUNKS=`cat /proc/cpuinfo |grep ^processor|wc -l`
	if [ "x$NCHUNKS" = "x" ]; then
	    NCHUNKS=1
	fi
    fi
else
    LSB_HOSTS_ARRAY=( $LSB_HOSTS )
    NCHUNKS=${#LSB_HOSTS_ARRAY[*]}
fi

# Magic to determine archiv method
case $in1 in
    *.gz)
	ZCAT=zcat
	;;
    *.bz2)
	ZCAT=bzcat
	;;
    *.xz)
	ZCAT=xzcat
	;;
    *)
	ZCAT=cat
esac
case $out1 in
    *.gz)
	GZIP=gzip
	;;
    *.bz2)
	GZIP=bzip2
	;;
    *.xz)
	GZIP=xz
	;;
    *)
	GZIP=cat
esac


#CUTADAPT_DIR=/data2/unige/dmg/snikolae/FB/cutadapt
#CUTADAPT=$CUTADAPT_DIR/bin/cutadapt
CUTADAPT=cutadapt


## Oh, do not dump to /tmp in multynode environment!
# tmpdir=/tmp/cutadapt.$USER.$$
tmpdir=cutadapt.$USER.$$
trap "rm -r $tmpdir" 0
mkdir $tmpdir


# Split input reads
$ZCAT $in1 > $tmpdir/in_1.fastq
READS=$[ `cat $tmpdir/in_1.fastq|wc -l` / 4 ]
SPLITLINES=$[ 4 * ( $READS / $NCHUNKS + $NCHUNKS ) ]

( split -l $SPLITLINES $tmpdir/in_1.fastq $tmpdir/in_1; rm $tmpdir/in_1.fastq )&
if [ -n "$in2" ]; then
    $ZCAT $in2 | split -l $SPLITLINES - $tmpdir/in_2 &
fi
wait


# Cut adapters
i=0
for f in $tmpdir/in_1*;do
    f2=$tmpdir/in_2${f#$tmpdir/in_1}
    of=$tmpdir/out_1${f#$tmpdir/in_1}
    of2=$tmpdir/out_2${f#$tmpdir/in_1}
    if [ -n "$in2" ]; then
	cmd="$CUTADAPT $* -f fastq -o $of -p $of2 $f $f2"
    else
	cmd="$CUTADAPT $* -f fastq -o $of $f"
    fi
    if [ "x$NO_LSB" = "xTrue" ]; then
	eval $cmd &
    else
	lsrun -v -m ${LSB_HOSTS_ARRAY[$i]} $cmd & 
    fi
    i=$[ $i + 1 ]
    if [ $i -gt $NCHUNKS ]; then
	echo "WHAT IS GOING ON?! More chunks than cores! $i>$NCHUNKS"
	i=0
    fi
done
wait

if [ "x$NO_LSB" != "xTrue" ]; then
    sleep 10 # Make sure NFS is in sync
fi

cat $tmpdir/out_1* | $GZIP >$out1 &
if [ -n "$in2" ]; then
    cat $tmpdir/out_2* | $GZIP >$out2 &
fi
wait
