import argparse
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from intervaltree import Interval, IntervalTree


parser = argparse.ArgumentParser( description='Parse FP')
parser.add_argument('-tab', '--tab', required=False)
parser.add_argument('-truvariFP', '--fp', required=False)

args = parser.parse_args()

FPList=[]
bcf_in = VariantFile(args.fp)

for rec in bcf_in.fetch():
    sv_id = rec.id
    FPList.append(rec.id)

tabFile = open(args.tab, "r")

for line in tabFile.readlines():
    id = line.split("_")[0]
    sv_len = line.split("_")[3]
    if id in FPList:
        print("%s\t%s" % (line.rstrip("\n"), "FP"))
    elif int(sv_len) >= 50:
        print("%s\t%s" % (line.rstrip("\n"), "TP"))
    else:
        print("%s\t%s" % (line.rstrip("\n"), "not_processed"))

