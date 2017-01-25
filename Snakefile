configfile: "config.json"

import os
TRANSCRIPTOME_INDEX = os.path.splitext(config["GENOME_GFF"])[0]+".tr"
SCRIPTDIR = os.path.join(os.path.dirname(os.path.abspath(__name__)), "scripts")


## Getting all samples in input directory and extensions
SAMPLES,SAMPLESEXT = glob_wildcards("input/{sample}.fastq{cext,|.(gz|bz2|xz)}")
if len(SAMPLES)==0:
    SAMPLES,steps = glob_wildcards("bams/{sample,[^.]*}.{steps}.bam")
    SAMPLES=list(set(SAMPLES))


def get_fastq_ext(sample):
    try:
        return SAMPLESEXT[SAMPLES.index(sample)]
    except:
        return ""


## Global targets
rule csv:
    input: expand("wigs/{sample}.stripped.trimmed.nonribo.csv", sample=SAMPLES)
rule bams:
    input: expand("bams/{sample}.stripped.trimmed.nonribo.bam", sample=SAMPLES)
rule randomcsv:
    input: expand("wigs/{sample}.stripped.trimmed.nonribo.random.csv", sample=SAMPLES)
rule randombams:
    input: expand("bams/{sample}.stripped.trimmed.nonribo.random.bam", sample=SAMPLES)
rule bamsnf:
    input: expand("bams/{sample}.stripped.trimmed.bam", sample=SAMPLES)
rule randombamsnf:
    input: expand("bams/{sample}.stripped.trimmed.random.bam", sample=SAMPLES)
rule csvnf:
    input: expand("wigs/{sample}.stripped.trimmed.csv", sample=SAMPLES)
rule randomcsvnf:
    input: expand("wigs/{sample}.stripped.trimmed.random.csv", sample=SAMPLES)


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


rule filter_ribosomal:
    input:
        fastq="cleaned/{sample}.fastq.gz",
        indexnc=config["INDEX_NC"]+".1.ebwt"
    output:
        "cleaned/{sample}.nonribo.fastq.gz"
    threads: 8
    shell:
    	"zcat {input.fastq} |"
	" bowtie -q -S -p $(( {threads}>1?{threads}-1:1 )) -v3 --seedlen=23 {config[INDEX_NC]} - |"
	" samtools  view -h -f4 - | samtools bam2fq - | gzip > {output}"


## I forgot, why I used fr-secondstrand
rule tophat:
    input:
        fastq="cleaned/{sample}.fastq.gz",
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

rule randomize_double_hits:
    input:
        "{sample}.bam"
    output:
        "{sample}.random.bam"
    shell:
        "{SCRIPTDIR}/randomizeDoubleHits.sh {input} {output}"

rule genWig:
    input:
        bam="bams/{sample}.bam"
    output:
        csv="wigs/{sample}.csv",
        pluscbg="wigs/{sample}.count.plus.bw",
        minuscbg="wigs/{sample}.count.minus.bw",
        # pluscw="wigs/{sample}.count.plus.wig",
        # minuscw="wigs/{sample}.count.minus.wig",
        plusbg="wigs/{sample}.plus.bw",
        minusbg="wigs/{sample}.minus.bw",
        # plusw="wigs/{sample}.plus.wig",
        # minusw="wigs/{sample}.minus.wig"
    script:
        "scripts/BamToBedGraph.R"

rule genRandomizedWig:
    input:
        bam="bams/{sample}.bam"
    output:
        csv="wigs/{sample,.*\.random}.csv",
        pluscbg="wigs/{sample}.count.plus.bw",
        minuscbg="wigs/{sample}.count.minus.bw",
        # pluscw="wigs/{sample}.count.plus.wig",
        # minuscw="wigs/{sample}.count.minus.wig",
        plusbg="wigs/{sample}.plus.bw",
        minusbg="wigs/{sample}.minus.bw",
        # plusw="wigs/{sample}.plus.wig",
        # minusw="wigs/{sample}.minus.wig"
    params:
        multihit=True
    script:
        "scripts/BamToBedGraph.R"

ruleorder: genRandomizedWig > genWig
