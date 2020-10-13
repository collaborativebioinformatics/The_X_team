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
    id = seqRecord.id.split("_")[0]
    seqLen = len(seqRecord.seq)
    LensDict[id] = seqLen


SVInfo = {}
AlLenHap1 = {}
AlLenHap2 = {}

for record in bam1:
    if record.is_secondary:
        continue
    if record.is_supplementary:
        continue
    if record.is_supplementary and record.is_secondary:
        continue
    id = record.query_name.split("_")[0]
    sv_len = record.query_name.split("_")[3]
    gt = record.query_name.split("_")[5]

    SVInfo[id] = [sv_len,gt]

    if int(record.flag) == 4:
        alignLen = 0
        AlLenHap1[id] = alignLen
    else:
        alignLen = int(record.reference_end) - int(record.reference_start)
        AlLenHap1[id] = alignLen

for rec in bam2:
    if rec.is_secondary:
        continue
    if rec.is_supplementary:
        continue
    if rec.is_supplementary and record.is_secondary:
        continue
    id2 = rec.query_name.split("_")[0]

    if int(rec.flag) == 4:
        alignLen2 = 0
        AlLenHap2[id2] = alignLen2
    else:
        alignLen2 = int(rec.reference_end) - int(rec.reference_start)
        AlLenHap2[id2] = alignLen2

for id in AlLenHap1.keys():
    if args.listSVs is not None and id in listSV:
        if id in FPList:
            print(id, LensDict[id], SVInfo[id][0], SVInfo[id][1], AlLenHap1[id], AlLenHap2[id], "FP")
        else:
            print(id, LensDict[id], SVInfo[id][0], SVInfo[id][1], AlLenHap1[id], AlLenHap2[id],"TP")

    if args.listSVs is None and int(SVInfo[id][0]) >= 50:
        if id in FPList:
            print(id, LensDict[id], SVInfo[id][0], SVInfo[id][1], AlLenHap1[id], AlLenHap2[id],"FP")
        else:
            print(id, LensDict[id], SVInfo[id][0], SVInfo[id][1], AlLenHap1[id], AlLenHap2[id],"TP")

    if args.listSVs is None and int(SVInfo[id][0]) < 50:
        if id in FPList:
            print(id, LensDict[id], SVInfo[id][0], SVInfo[id][1], AlLenHap1[id], AlLenHap2[id],"FP-le50")
        else:
            print(id, LensDict[id], SVInfo[id][0], SVInfo[id][1], AlLenHap1[id], AlLenHap2[id],"le50")


#for b in bam.fetch(chrom, startLook, endLook, until_eof=True):