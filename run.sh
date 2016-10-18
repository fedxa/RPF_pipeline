#!/bin/sh
dir=`dirname $0`

# module add Utility/snakemake/3.5.4
# Note -- do not use snamakemake 3.8.0 -- race condition.  3.7.x is ok (install pyaml also)

CLUSTER=cluster.json


module add UHTS/Aligner/bowtie/0.12.9
module add UHTS/Aligner/tophat/2.0.13
module add UHTS/Analysis/samtools/1.3
module add UHTS/Quality_control/fastqc/0.11.2
module add UHTS/Analysis/fastx_toolkit/0.0.13.2
module add R/3.2.2


snakemake --cluster-config $CLUSTER \
	  --cluster-sync "bsub -K -q {cluster.queue} -o {output[0]}.log -n {threads} -M {cluster.mem} -R 'span[hosts=1]'" \
	  --jobs ${JOBS:=100} \
	  -p $*
