#!/bin/bash

## Wrapper script to do the merging block in the diploid pipeline workflow to merge assembly errors due to undercollapsed heterozygosity.
## Requires LASTZ aligner in path
## Written by: OD

USAGE="
*********************************
 bash run-asm-merger.sh <cprops> <asm> <faSplitFolder>
*********************************
"

## HANDLE OPTIONS
while getopts "b:h" opt; do
case $opt in
	b)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
            echo ":) -b flag was triggered, assuming ${OPTARG} wide band alternative haplotype detection" >&1
            band_size=$OPTARG
        else
            echo ":( Wrong syntax for band size. Using the default value ${band_size}" >&2
        fi
	;;
	h) echo "$USAGE"
		exit 0
	;;
	*) echo "$USAGE" >&2
		exit 1
	;;
esac
done
shift $(( OPTIND-1 ))


# Handle dependencies
pipeline=`cd "$( dirname $0)" && cd .. && pwd`

# Handle arguments
cprops=$1
asm=$2
faSplit=$3

#	1) Perform pariwise alignment in a band of predefined size (default: 1Mb) and output reliable overlap input and original_asm annotation track (overlaps_2D_input.txt and overlaps_2D_asm.txt)
	bash ${pipeline}/merge/run-pairwise-alignment.sh -b ${band_size} ${cprops} ${asm} ${faSplit}


#	2) Extract connected components, vote orientation, tile cluster, vote order
	awk -f ${pipeline}/merge/identify-clusters-w-orientation.awk overlaps_2D_input.txt | awk -f ${pipeline}/merge/vote-orientation.awk ${cprops} ${asm} - | awk -f ${pipeline}/merge/tile-cluster-insides.awk ${cprops} overlaps_2D_input.txt - | awk -f ${pipeline}/merge/vote-order.awk ${cprops} - ${asm} > `basename ${asm} .asm`"_tiled.asm"

	tiled_asm=`basename ${asm} .asm`"_tiled.asm"
	
#	3) Merge tiled asm
	bash ${pipeline}/merge/merge-tiled-asm.sh ${cprops} ${tiled_asm} ${faSplit}

	

