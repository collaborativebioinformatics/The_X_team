#!/usr/bin/env bash

sort -k4,4 $1 > $1.sorted
sort -k4,4 $2 > $2.sorted

join -1 4 -2 4 $1.sorted $2.sorted | tr " " "\t" | grep -v chrom > $3 
