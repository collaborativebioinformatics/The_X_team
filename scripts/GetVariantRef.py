import argparse
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np


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
    sv_end = rec.stop
    sv_reads = rec.info['RNAMES']
    sv_read_support = rec.info['RE']
    sv_len = rec.info['SVLEN']
    sv_type = rec.info['SVTYPE']
    sv_filter = rec.filter
    print(sv_id,sv_chr,sv_start,sv_end)

