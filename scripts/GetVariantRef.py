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

    sv_type = rec.info['SVTYPE']

    if "INS" in rec.id or "DEL" in rec.id:
        sv_id = rec.id
    else:
        sv_id = str(rec.id) + "." + str(sv_type)

    sv_chr = rec.chrom
    sv_start = rec.pos
    sv_len = abs(rec.info['SVLEN'])
    sv_type = rec.info['SVTYPE']
    sv_read_support = rec.info['RE']

    gt_info_cuteSV = str(rec.samples[0]['GT'][0]) + "/" + str(rec.samples[0]['GT'][1])

    gt_info = ""
    if gt_info_cuteSV == "0/0":
        gt_info = "Lowfreq"
    elif gt_info_cuteSV == "0/1":
        gt_info = "Het"
    elif gt_info_cuteSV == "1/1":
        gt_info = "Hom"
    else:
        gt_info = "None"

    for fil in rec.filter:
        sv_filter = fil

    sv_end = int(sv_start) + int(sv_len)
    seq_ID = str(sv_id) + "_" + str(sv_chr) + "_" + str(sv_start) + "_" + str(sv_end)

    dictSVcoords[sv_id] = [sv_chr, sv_start, sv_end, sv_filter, gt_info, sv_read_support]

    sv_allele = rec.alleles[1]

    #print(sv_id, sv_chr, sv_start, sv_end, rec.alleles)

    if "INS" in sv_type:
        dictAllelesINS[sv_id] = sv_allele

dicRef = refParser(args.reference)

treeSV = clusterNearSVs(dictSVcoords)


for SV in dictSVcoords.keys():
    coords = dictSVcoords[SV]

    sv_chr = coords[0]
    sv_start = coords[1]
    sv_start_prev = int(coords[1]) - 1
    sv_start_next = int(coords[1]) + 1
    sv_end = coords[2]
    sv_len = sv_end - sv_start

    sv_filter = coords[3]
    gt_info = coords[4]
    sv_read_support = coords[5]

    ref_start_toExtract = sv_start - 1500
    ref_end_toExtract = sv_end + 1500

    if "INS" in SV:
        InsertedSeq = dictAllelesINS[SV]
        print(">%s_%s_%s_%s_%s_%s_%s_ref%s-%s\n%s%s%s" % (SV, sv_chr, sv_start, sv_len, sv_filter, gt_info, sv_read_support, ref_start_toExtract, ref_end_toExtract, dicRef[sv_chr][ref_start_toExtract:sv_start_prev], dictAllelesINS[SV], dicRef[sv_chr][sv_start_next:ref_end_toExtract]))
    elif "DEL" in SV:
        DeletedSeq = dicRef[sv_chr][sv_start:sv_start+sv_len]
        print(">%s_%s_%s_%s_%s_%s_%s_ref%s-%s\n%s%s" % (SV, sv_chr, sv_start, sv_len, sv_filter, gt_info, sv_read_support, ref_start_toExtract,ref_end_toExtract, dicRef[sv_chr][ref_start_toExtract:sv_start_prev],dicRef[sv_chr][sv_start+sv_len:ref_end_toExtract]))