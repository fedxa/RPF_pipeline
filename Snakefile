configfile: "config.json"

import os
TRANSCRIPTOME_INDEX = os.path.splitext(config["GENOME_GFF"])[0]+".tr"
SCRIPTDIR = os.path.join(os.path.dirname(os.path.abspath(__name__)), "scripts")


## Getting all samples in input directory and extensions
SAMPLES,SAMPLESEXT = glob_wildcards("input/{sample}.fastq{cext,|.(gz|bz2|xz)}")

def get_fastq_ext(sample):
    try:
        return SAMPLESEXT[SAMPLES.index(sample)]
    except:
        return ""


## Global targets
rule bams:
    input: expand("bams/{sample}.bam", sample=SAMPLES)
rule csv:
    input: expand("wigs/{sample}.csv", sample=SAMPLES)


## Utility rules: generating indices
rule bowtie_index:
    input:
        "{index}.fa"
    output:
        "{index}.1.ebwt"
    shell:
        "bowtie-build {input} {wildcards.index}"

rule tophat_transcriptome_index:
    input:
        gff="{genome}.gff",
        index=config["INDEX"]+".1.ebwt"
    output:
        "{genome}.tr.1.ebwt"
    shell:
        "tophat --bowtie1 -G {input.gff} --transcriptome-index={wildcards.genome}.tr {config[INDEX]}"

ruleorder: tophat_transcriptome_index > bowtie_index

rule faidx:
    input:
        "{index}.fa"
    output:
        "{index}.fa.fai"
    shell:
        "samtools faidx {input}"


## Main pipeline
rule cutadapt:
    input:
     	lambda wildcards: "input/{sample}.fastq{ext}".format(sample=wildcards.sample, ext=get_fastq_ext(wildcards.sample))
    output:
        "cleaned/{sample}.stripped.fastq.gz"
    threads: 16
    shell:
        "NCHUNKS={threads} {SCRIPTDIR}/cutadapt_parallel.sh"
        " {input} '' {output} '' -a {config[ADAPTER]} --trimmed-only -m 21 -e 0.1"

rule trim_bad_ligation:
    input:
        "cleaned/{sample}.stripped.fastq.gz"
    output:
        "cleaned/{sample}.stripped.trimmed.fastq.gz"
    threads: 2
    shell:
        "zcat {input} | fastx_trimmer -Q33 -f 2 | gzip > {output}"

def sample_trimmed_or_not(wildcards):
    "Return the trimmed or untrimmed FASTQ file name depending in TRIM_BAD_LIGATION config"
    if config["TRIM_BAD_LIGATION"]:
        return "cleaned/{sample}.stripped.trimmed.fastq.gz".format(sample=wildcards.sample)
    else:
        return "cleaned/{sample}.stripped.fastq.gz".format(sample=wildcards.sample)

## I forgot, why I used fr-secondstrand
rule tophat_unfiltered:
    input:
        fastq=sample_trimmed_or_not,
        index=config["INDEX"]+".1.ebwt",
        index_tr=TRANSCRIPTOME_INDEX+".1.ebwt"
    output:
        dir="bams/{sample}",
        accepted="bams/{sample}/accepted_hits.bam"
    threads: 16
    version: 0.1
    shell:
        "tophat -p{threads} --bowtie1"
	" --transcriptome-index={TRANSCRIPTOME_INDEX} --no-novel-juncs"
	" --library-type fr-secondstrand"
        " -o {output.dir} {config[INDEX]} {input.fastq}"

## Fix them form the strand in the alignment flag
## https://genomebytes.wordpress.com/2013/07/08/fixing-the-xs-tag-in-tophat-output-bug-fixing/        
rule fix_tophat_output:
    input:
        "{sample}/accepted_hits.bam"
    output:
        "{sample}.bam"
    threads: 3
    shell:
        "{SCRIPTDIR}/fix_tophat_XS_tags.sh {input} > {output}"

rule genWig:
    input:
        bam="bams/{sample}.bam",
        config="config.json",
    output:
        csv="wigs/{sample}.csv",
        plusbg="wigs/{sample}.plus.bedgraph.gz",
        minusbg="wigs/{sample}.minus.bedgraph.gz"
    script:
        "scripts/BamToBedGraph.R"

rule bedGraphToBigWig:
    input:
        bg="{file}.bedgraph.gz",
        fai=config["INDEX"]+".fa.fai"
    output:
        "{file}.bw"
    shell:
        "bedGraphToBigWig {input.bg} {input.fai} {output}"
