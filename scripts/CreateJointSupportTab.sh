#!/usr/bin/env bash

paste <(echo "ID")  <( head -1 $1) <(head -1 $2 | cut -f 2- | awk '{ s=""; for (i=1;i<NF;i++) { s=s "2_"$i ; if (i<NF-1) { s=s" "; } } print s }') > $3
tail -n +2 $1 | sort -k4,4  > $1.sorted
tail -n +2 $2 | sort -k4,4  > $2.sorted

join -1 4 -2 4 $1.sorted $2.sorted | tr " " "\t"  >> $3
