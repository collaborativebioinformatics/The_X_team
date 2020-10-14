import argparse
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from intervaltree import Interval, IntervalTree

def readBAM(bam):
    ChrDict = {}
    for record in bam:
        id = record.query_name.split("_")[0]
        start = record.reference_end
        chr = record.reference_name
        ChrDict[id] = [chr,start]
    return ChrDict


parser = argparse.ArgumentParser( description='Parse FP')
parser.add_argument('-tab', '--tab', required=False)
parser.add_argument('-truvariFP', '--fp', required=False)
parser.add_argument('-truvariFP-Param', '--fpParam', required=False)
parser.add_argument('-hap1', '--hap1', required=False)
parser.add_argument('-hap2', '--hap2', required=False)

args = parser.parse_args()

FPList=[]
bcf_in = VariantFile(args.fp)

FPListExtra=[]
bcf_in_param = VariantFile(args.fpParam)

for rec in bcf_in.fetch():
    sv_id = rec.id
    FPList.append(rec.id)

for rec in bcf_in_param.fetch():
    sv_id = rec.id
    FPListExtra.append(rec.id)

bam1 = pysam.AlignmentFile(args.hap1, "rb", header=True)
bam2 = pysam.AlignmentFile(args.hap2, "rb", header=True)

ChrHap1 = readBAM(bam1)
ChrHap2 = readBAM(bam2)

tabFile = open(args.tab, "r")

for line in tabFile.readlines():
    id = line.split("_")[0]
    sv_len = line.split("_")[3]
    if id in FPList and id in FPListExtra:
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (line.rstrip("\n"), "FP", "FP-Param", ChrHap1[id][0],  ChrHap1[id][1], ChrHap2[id][0],  ChrHap2[id][1] ))
    if id in FPList and id not in FPListExtra:
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (line.rstrip("\n"), "FP", "TP-Param", ChrHap1[id][0],  ChrHap1[id][1], ChrHap2[id][0],  ChrHap2[id][1] ))
    elif int(sv_len) >= 50:
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (line.rstrip("\n"), "TP", "TP-Param", ChrHap1[id][0],  ChrHap1[id][1], ChrHap2[id][0],  ChrHap2[id][1] ))
    else:
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (line.rstrip("\n"), "not_processed", "not_processed", ChrHap1[id][0],  ChrHap1[id][1], ChrHap2[id][0],  ChrHap2[id][1] ))

