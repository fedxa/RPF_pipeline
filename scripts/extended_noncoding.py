#!/usr/bin/env python3
#
# For all features in the index of the noncoding RNA fasta
# Get the fearues extended by _2_ nucleotide beyond start (end for - strand features)
#
# Check, that that can be doen by selecting all noncoding_exons from the genome gff?

import pandas as pd
import io
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.Alphabet import generic_dna
import argparse


def read_genome(fasta, full=True, annotate=False):
    """Read the genome, adding annotations if suggested.

    file -- gff file to read
    full -- read the whole genome into memory (default=True)
    annotate -- read the annotations (forces full=True)"""
    if full or annotate:
        with open(fasta, 'r') as fd:
            yg = SeqIO.to_dict(SeqIO.parse(fd, "fasta", generic_dna))
        if annotate:
            create_annotations(yg, annotate)
    else:
        yg = SeqIO.index(fasta, "fasta")
    return yg


def get_field_from_fields(t, field):
    res = t.extra_fields.str.extract('(?:^|;)'+field+'=([^;]*)')
    return res.str.replace('%20', ' ')


def fields2dict(id):
    return {x[0]: x[1] for x in [r.split('=') for r in id.split(';')]}


def get_all_keys_in_fields(ann):
    kk = set()
    for k in ann.extra_fields.map(lambda x: fields2dict(x).keys()):
        for el in k:
            kk.add(el)
    return kk


def read_gff_annotations(file):
    """Reads pandas DataFrame from GFF file, extending the attributes into
    separate columns"""
    buffer = ''
    with open(file, 'r') as f:
        for l in f:
            if l.startswith('###'):    # Start of FASTA part, stop
                break
            if not l.startswith('#'):  # Drop comments
                buffer += l
    buffd = io.StringIO(buffer)
    res = pd.read_csv(buffd, sep='\t',
                      names=['chr', 'source', 'type', 'start',
                             'end', 'score', 'strand', 'frame',
                             'extra_fields'],
                      na_values={'source': '.', 'score': '.', 'strand': '?'})
    for f in get_all_keys_in_fields(res):
        res[f] = get_field_from_fields(res, f)
    del(res['extra_fields'])
    return res


def read_regions(file):
    """Read names of the sequences in FASTA file"""
    res = []
    with open(file, 'r') as f:
        for l in f:
            if l.startswith('>'):
                wd = l.split()
                res += [ wd[0][1:] ]
    return res


def select(gff, selected):
    return gff[gff.ID.isin(selected[0])]


def find_feature(name, ann):
    res = ann[ann['ID'] == name]
    if (res.empty):
        return None
    else:
        if (len(res) > 1):
            print('Warning: multiple genes named "'+name+'"!')
        return res.iloc[0]


def gene_seq(gene, ann, yg, part='CDS'):
    loc = gene_feature(gene, ann, part)
    return loc.extract(yg[gene['chr']])


def gene_feature(gene, ann, part):
    """Get the compound location for the CDS for the gene"""
    if ((part == 'all')):
        return FeatureLocation(gene['start']-1, gene['end'])
    else:
        cds = ann[(ann.Parent == gene['ID']) & (ann.type == part)]
        return cds.apply(lambda e: FeatureLocation(e['start']-1, e['end']),
                         axis=1).sum()


def get_extended_feature(f, ann, extendby=3, extendby_end=3):
    gf = gene_feature(f, ann, 'noncoding_exon')
    if f.strand == '-':
        gf = FeatureLocation(int(f['start']-1-extendby_end), int(f['start']-1))+gf+FeatureLocation(int(f['end']), int(f['end']+extendby))
    else:
        gf = FeatureLocation(int(f['start']-1-extendby), int(f['start']-1))+gf+FeatureLocation(int(f['end']), int(f['end']+extendby_end))
    return gf


def get_all_features(ann, selected, yg, extendby=3, extendby_end=3):
    sequences = []
    for fn in selected:
        ef = find_feature(fn, ann)
        gf = get_extended_feature(ef, ann, extendby=extendby, extendby_end=extendby_end)
        seq = gf.extract(yg[ef.chr])
        if ef.strand == '-':
            seq = seq.reverse_complement()
        seq.name = fn
        seq.id = fn
        sequences.append(seq)
    return sequences


def save_sequences(seqs, fname):
    with open(fname, 'w') as fd:
        SeqIO.write(seqs, fd, 'fasta')

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extend noncoding regions")
    parser.add_argument('-b', '--begin', type=int, default=3,
                        help='Extend from beginning by')
    parser.add_argument('-e', '--end', type=int, default=3,
                        help='Extend from end by')
    parser.add_argument('-o', '--outfile',
                        default='data/Yeast-Noncoding-extended.fa',
                        help="Output FASTA file")
    parser.add_argument('infile', nargs='?',
                        default='data/Yeast-Noncoding.fa',
                        help="SAM file")
    args = parser.parse_args()

    ann = read_gff_annotations('data/yeast.gff')
    yg = read_genome('data/yeast.fa')
    sel = read_regions(args.infile)
    selseqs = get_all_features(ann, sel, yg,
                               extendby=args.begin, extendby_end=args.end)
    save_sequences(selseqs, args.outfile)
