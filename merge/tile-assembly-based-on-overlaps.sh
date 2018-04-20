#!/bin/bash
## Sandboxing part of the merge pipeline that performs connected component analysis and tiling based on pairwise alignment data
## USAGE: bash ./merge/tile-assembly-based-on-overlaps.sh <cprops> <asm> <overlaps-from-alignment-2D-annotation-file>
## Input: cprops, asm, overlaps 2D annotation files
## Output: _tiled.asm file, an annotated .asm file in which connected components are annotated with {...} 
## Steps involved: extract connected components, vote orientation, tile cluster, vote order
## Written by: OD
## Version: 180402

# Handle dependencies
pipeline=`cd "$( dirname $0)" && cd .. && pwd`

# Handle arguments
cprops=$1
asm=$2
overlaps=$3

([ -z $1 ] || [ -z $2 ] || [ -z $3 ]) && echo >&2 "Not sure how to parse your input: files not listed or not found at expected locations. Exiting!" && echo >&2 "$USAGE_short" && exit 1

awk -f ${pipeline}/merge/identify-clusters-w-orientation.awk ${overlaps} | awk -f ${pipeline}/merge/vote-orientation.awk ${cprops} ${asm} - | awk -f ${pipeline}/merge/tile-cluster-insides.awk ${cprops} ${overlaps} - | awk -f ${pipeline}/merge/vote-order.awk ${cprops} - ${asm} > `basename ${asm} .asm`"_tiled.asm"
