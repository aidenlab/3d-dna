#!/bin/bash

## Wrapper script to do the merging block in the diploid pipeline workflow to merge assembly errors due to undercollapsed heterozygosity.
## Requires LASTZ aligner in path
## The idea is that sequence similarities identified via LASTZ alignments need to be long enough and with high enough identity score to meet a certain (length x identity) threshold. An additional (very high) alignment score threshold is listed that represents saturation of the length x identity condition.
## Written by: OD
## Version: 170217

USAGE="
*********************************
 bash run-asm-merger.sh [options] <cprops> <asm> <faSplitFolder>
*********************************
"

# defaults
merger_search_band=3000000
merger_alignment_score=50000000
merger_alignment_identity=20
merger_alignment_length=20000
merger_lastz_options=\"--gfextend\ --gapped\ --chain=200,200\"


## HANDLE OPTIONS
while getopts "b:s:i:l:o:t:h" opt; do
case $opt in
	b)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
            echo ":) -b flag was triggered, assuming ${OPTARG} bp-wide band for alternative haplotype detection." >&1
            merger_search_band=$OPTARG
        else
            echo ":( Wrong syntax for band size. Using the default value ${merger_search_band}." >&2
        fi
	;;
	l) re='^[0-9]+$'
		if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
            echo ":) -l flag was triggered, assuming ${OPTARG} bp as minimal overlap length at which sequences might be classified as alternative haplotypes." >&1
            merger_alignment_length=$OPTARG
        else
            echo ":( Wrong syntax for acceptable synteny length. Using the default value ${merger_alignment_length}." >&2
        fi
	;;
    i) re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
            echo ":) -i flag was triggered, assuming ${OPTARG} as minimal identity (per length) at which similar sequences might be classified as alternative haplotypes." >&1
            merger_alignment_identity=$OPTARG
        else
            echo ":( Wrong syntax for acceptable identity score. Using the default value ${merger_alignment_identity}." >&2
        fi
    ;;
	s) re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
            echo ":) -s flag was triggered, assuming ${OPTARG} as saturated alignment score at which similar sequences are classified as alternative haplotypes, irrespective of length." >&1
            merger_alignment_score=$OPTARG
        else
            echo ":( Wrong syntax for acceptable alignment score. Using the default value ${merger_alignment_score}." >&2
        fi
    ;;
	o) re='^\"--.+\"$'
		if [[ $OPTARG =~ $re ]]; then
			echo ":) -o flag was triggered, assuming ${OPTARG} as a list of options to pass on to LASTZ to tune alignment." >&1
			merger_lastz_options=${OPTARG}
		else
			echo ":( Wrong syntax for LASTZ option string. Using the default value ${merger_lastz_options}." >&2
		fi
	;;
	t)
		if [ -f ${OPTARG} ]; then
			echo ":) -t flag was triggered, assuming $OPTARG contains tiling results and fast-forwarding to merging per se." >&1
			tiled_asm=$OPTARG
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

if [ -z ${tiled_asm} ]; then
#	1) Perform pariwise alignment in a band of predefined size (default: 1Mb) and output reliable overlap input and original_asm annotation track (overlaps_2D_input.txt and overlaps_2D_asm.txt)
	bash ${pipeline}/merge/align-nearby-sequences-and-filter-overlaps.sh -b ${merger_search_band} -s ${merger_alignment_score} -i ${merger_alignment_identity} -l ${merger_alignment_length} -o "${merger_lastz_options}" ${cprops} ${asm} ${faSplit}

#	2) Extract connected components, vote orientation, tile cluster, vote order
	bash ${pipeline}/merge/tile-assembly-based-on-overlaps.sh ${cprops} ${asm} overlaps_2D_input.txt
	tiled_asm=`basename ${asm} .asm`"_tiled.asm"
fi

#	3) Merge tiled asm
	bash ${pipeline}/merge/merge-tiled-asm.sh -o "${merger_lastz_options}" ${cprops} ${tiled_asm} ${faSplit}

	

