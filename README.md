Ribosome footprinting pipeline
==============================

This is a simple [Snakemake][] based implementation of a Ribosome
profiling pipeline


Installation
------------

### Requirements

- R > 3.3.1
  + GenomicAlignments
  + GenomicFeatures
  + RiboProfiling
  + rtracklayer
- tophat
- bowtie
- FastX
- samtools

Run
---


Configuration
-------------

The configuration is in the file `config.json` which looks like

    {
    "INDEX": "data/yeast",
    "GENOME_GFF": "data/yeast.gff",
    "ADAPTER": "CTGTAGGCACCATCAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    "MANUAL_SHIFTS_FILE": "shifts.csv"
    }

The entries mean:

  * *INDEX* basename of the FASTA file with the reference sequence.
    In the example above the reference should be named `data/yeast.fa`
  * *GENOME_GFF* name of the genome annotation in gff format
  * *ADAPTER* Adapter to strip of the end of the sequence
  * *MANUAL_SHIFTS_FILE* Name of the csv with the manual selection of
    the reads and position of the P-sites. See P-site WIGs
	

P-site WIG generation
---------------------

The WIG (or, currently, BedGraph) files are being generated by
selecting only the reads with length that is close to the correct one,
and then selecting the fixed offset into the read, which is individual
for each read length. By default, the algorithm selects the length
corresponding to the maximum of the length histogram +- one.  The
histogram with color-marked selected lengths is saved in
`wig/{sample}_Length.pdf`.  Then, for each of the lenght individually
the metagene distribution of the starts of the reads relative to the
start of the nearest OFR is plotted in `wig/{sample}_Asiteshifts.pdf`.
(_TODO:_ is it A or P sites? recheck the terminology.) Then, the
corrsponding offset into the read is taken, which leads to the
resulting bedgraph files.

However, if this fails, it is possible to create the
*MANUAL_SHIFTS_FILE* with the sete of the read lengths and shifts for
each sample. It has three columns: file name, length list (space
separated), and shift list. the pipeline creates the CSV files with
the values it found automatically in `wigs/{sample}.shifts.csv`, which
can be used as the template to create the *MANUAL_SHIFTS_FILE*.


Notes on linkers/adapters in Ribosome Footprinting
--------------------------------------------------

### Detailed sequence breakout for Ribosome footprint library

                                                                                          / This one (rApp) seems to not be present!
    AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTATGCATGCATGCATGCATGCATGCATGCaCTGTAGGCACCATCAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG
                                                              ^^^^Sequence^^^^^^^^^^^^^^^^-CTGTAGGCACCATCAATAGATCGGAAGAGC                     TGACCA
    AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT                                               GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
    TruSeq Universal Adapter                                                                                 TruSeq Indexed Adapter           Index
    AATGATACGGCGACCACCGA
    Portion bind to flow cell
                                                                                                                                                        CGTATGCCGTCTTCTGCTTG
                                                                                                                                                        Binds to flow cell for paired end reads
                                                                                                             GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -> Index read proceeds
                                                                                                             Index read primer
                             ACACTCTTTCCCTACACGACGCTCTTCCGATCT -> Read proceeds
                             Read 1 primer

### Linkers (some references)

1. [Adapter](http://onlinelibrary.wiley.com/doi/10.1002/wrna.1172/pdf)

   `5'rAppCTGTAGGCACCATCAAT/3ddC/3'`

2. [cloning linkers](https://eu.idtdna.com/pages/landing/cloning-linkers)

	`AATGATACGGCGACCACCGAGATCTACAC`  
	`CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCAGACGTGTGCTCTTCCG`: IDX4


[Snakemake]: https://bitbucket.org/snakemake/snakemake/wiki/Home
