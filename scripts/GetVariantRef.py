import argparse
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def refParser(ref):
    dictRef = {}
    fastaParse = SeqIO.parse(ref, "fasta")
    for seqRecord in fastaParse:
        chr = seqRecord.id
        seq = seqRecord.seq
        dictRef[chr] = seq
    return dictRef


parser = argparse.ArgumentParser( description='Parse FP')
parser.add_argument('-vcf', '--vcf', required=False)
parser.add_argument('-ref', '--reference', required=False)
parser.add_argument('-o', '--out', required=False)

args = parser.parse_args()

fname = args.vcf
bcf_in = VariantFile(fname)

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

    print(sv_id,sv_chr,sv_start, sv_start, sv_len)

    sv_end = int(sv_start) + int(sv_len)
    seq_ID = str(sv_id) + "_" + str(sv_chr) + "_" + str(sv_start) + "_" + str(sv_end)

    print(sv_end,seq_ID)


dicRef = refParser(args.reference)

for chr in dicRef.values():
    print(chr[1])