#!/usr/bin/env python

import pysam
import argparse 
import bisect

ap=argparse.ArgumentParser(description="Score alignments for support of SV")

ap.add_argument("--bam", help="Aligned contigs")
ap.add_argument("--window", help="padding window", type=int, default=1500)
ap.add_argument("--skip", help="Amount of bases to skip to avoid full pairwise alignment", type=int, default=20)

args=ap.parse_args()

bam=pysam.AlignmentFile(args.bam)

def ParseHeader(name):
    name=name.replace("low_freq", "LowFreq")

    vals=name.split("_")
    svtype=vals[0].split(".")[1]
    return {"name" : vals[0], 
            "chrom" : vals[1],
            "start" : vals[2],
            "svlen" : vals[3],
            "gt" : vals[4],
            "svtype" : svtype   }


print("ID\tqLen\ttSpan\tminDiff\tsvtype\tsvlen")
for aln in bam:
    rEnd   = aln.reference_end
    rStart = aln.reference_start
    if (rEnd - rStart < 2*args.window):
        continue
    pairs=aln.get_aligned_pairs(matches_only=True)
    qPos=[p[0] for p in pairs]
    rPos=[p[1] for p in pairs]

    qStart = bisect.bisect_left(qPos, args.window-20)
    qEnd = min(args.window+10,bisect.bisect_left(qPos, qPos[-1]-args.window+20))

    tooBig=100000000
    minDiff = tooBig
    minS = None
    minE = None
    for s in range(qStart, 0, -args.skip):
        for e in range(qEnd, len(qPos), args.skip):
            qLen=qPos[e] - qPos[s]
            rLen=rPos[e] - rPos[s]
            alnDiff = abs(qLen-rLen)
            if alnDiff < minDiff:
                minDiff = alnDiff
                minS=s
                minE=e
    header=ParseHeader(aln.query_name)

    if aln.seq is not None:
        print(aln.query_name + "\t" + str(len(aln.seq)) + "\t" + str(rEnd-rStart) + "\t" + str(minDiff) + "\t" + header["svtype"] + "\t" + header["svlen"])

    

