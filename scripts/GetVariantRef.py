import argparse
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from intervaltree import Interval, IntervalTree

def refParser(ref):
    dictRef = {}
    fastaParse = SeqIO.parse(ref, "fasta")
    for seqRecord in fastaParse:
        chr = seqRecord.id
        seq = seqRecord.seq
        dictRef[chr] = seq
    return dictRef


def clusterNearSVs(dictSVcoords):
    trees = {}
    for SV in dictSVcoords.keys():
        coords = dictSVcoords[SV]
        sv_chr = coords[0]
        sv_start = int(coords[1]) - 1500
        sv_end = int(coords[2]) + 1500
        if sv_chr not in trees:
            trees[sv_chr] = IntervalTree()
        trees[sv_chr][sv_start:sv_end] = str(SV)
    return trees


parser = argparse.ArgumentParser( description='Parse FP')
parser.add_argument('-vcf', '--vcf', required=False)
parser.add_argument('-ref', '--reference', required=False)
parser.add_argument('-o', '--out', required=False)

args = parser.parse_args()

fname = args.vcf
bcf_in = VariantFile(fname)

dictSVcoords = {}
dictAllelesINS = {}

for rec in bcf_in.fetch():
    if "BND" in str(rec):
        continue
    sv_id = rec.id
    sv_chr = rec.chrom
    sv_start = rec.pos
    sv_reads = rec.info['RNAMES']
    sv_read_support = rec.info['RE']
    sv_len = abs(rec.info['SVLEN'])
    sv_type = rec.info['SVTYPE']
    sv_filter = rec.filter

    sv_end = int(sv_start) + int(sv_len)
    seq_ID = str(sv_id) + "_" + str(sv_chr) + "_" + str(sv_start) + "_" + str(sv_end)

    dictSVcoords[sv_id] = [sv_chr, sv_start, sv_end]

    sv_allele = rec.alleles[1]

    if "INS" in sv_type:
        dictAllelesINS[sv_id] = sv_allele

dicRef = refParser(args.reference)

treeSV = clusterNearSVs(dictSVcoords)

print(dicRef["1"][10369:10373])

for SV in dictSVcoords.keys():
    coords = dictSVcoords[SV]

    sv_chr = coords[0]
    sv_start = coords[1]
    sv_start_prev = int(coords[1]) - 1
    sv_start_next = int(coords[1]) + 1
    sv_end = coords[2]
    sv_len = sv_end - sv_start

    ref_start_toExtract = sv_start - 1500
    ref_end_toExtract = sv_end  + 1500

    if "INS" in SV:
        print(SV, sv_chr, sv_start, sv_end, dicRef[sv_chr][ref_start_toExtract:sv_start_prev],dictAllelesINS[SV],dicRef[sv_chr][sv_start_next:ref_end_toExtract])
    elif "DEL" in SV:
        print(SV, sv_chr, sv_start, sv_end, dicRef[sv_chr][ref_start_toExtract:sv_start_prev],dicRef[sv_chr][sv_start_next+sv_len:ref_end_toExtract])