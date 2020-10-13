import argparse
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from intervaltree import Interval, IntervalTree


parser = argparse.ArgumentParser( description='Parse FP')
parser.add_argument('-hap1', '--hap1', required=False)
parser.add_argument('-hap2', '--hap2', required=False)
parser.add_argument('-varfasta', '--varfasta', required=False)
parser.add_argument('-vcf', '--vcf', required=False)
parser.add_argument('-truvariFP', '--fp', required=False)
parser.add_argument('-o', '--out', required=False)
parser.add_argument('-list', '--listSVs', required=False)
args = parser.parse_args()


bam1 = pysam.AlignmentFile(args.hap1, "rb", header=True)
bam2 = pysam.AlignmentFile(args.hap2, "rb", header=True)

fastaParse = SeqIO.parse(args.varfasta, "fasta")

if args.listSVs is not None:
    filSV = open(args.listSVs, "r")
    listSV=[]
    for line in filSV.readlines():
        listSV.append(line.rstrip("\n"))

FPList=[]
bcf_in = VariantFile(args.fp)

for rec in bcf_in.fetch():
    sv_id = rec.id
    FPList.append(rec.id)


LensDict = {}

for seqRecord in fastaParse:
    id = seqRecord.id
    seqLen = len(seqRecord.seq)
    LensDict[id] = seqLen


for record in bam1:
    id = record.query_name.split("_")[0]
    sv_len = record.query_name.split("_")[3]

    if int(record.flag) != 4:
        alignLen = int(record.reference_end) - int(record.reference_start)
        if record.is_secondary:
            continue
        if record.is_supplementary:
            continue
        if record.is_supplementary and record.is_secondary:
            continue
        if args.listSVs is not None and id in listSV:
            if id in FPList:
                print(record.query_name, LensDict[record.query_name], alignLen,"FP")
            else:
                print(record.query_name, LensDict[record.query_name], alignLen,"TP")

        if args.listSVs is None and int(sv_len) >= 50:
            if id in FPList:
                print(record.query_name, LensDict[record.query_name], alignLen,"FP")
            else:
                print(record.query_name, LensDict[record.query_name], alignLen,"TP")

        if args.listSVs is None and int(sv_len) < 50:
            if id in FPList:
                print(record.query_name, LensDict[record.query_name], alignLen,"FP-le50")
            else:
                print(record.query_name, LensDict[record.query_name], alignLen,"le50")

#for b in bam.fetch(chrom, startLook, endLook, until_eof=True):