#!/bin/bash
## Wrapper script to analyze the assembly and place fragments from among the small scaffolds back into the assembly
## NOTE: now is run after splitting but for diploid pipeline could be run for for tiled assembly
## NOTE: Relies on standard annotations :::fragment_ and :::debris!
## NOTE: Probably should be done after chrom splitting to avoid a slim chance that a misassembled contig/scaffold spanned chromosomes, and the chromosomes ended up in the order and orientation to fit the edges of the contig/scaffold

USAGE="
**********************************
	
	./seal-asm.sh <current_cprops> <current_asm>
	
**********************************
"

## HANDLE OPTIONS

while getopts "s:h" opt; do
case $opt in
	h) echo "$USAGE"
		exit 0
	;;
	s) 	re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
	          echo "... -s flag was triggered, will attempt to place back only singleton debris contigs/scaffolds and those less than $OPTARG" >&1
	          SIZE=$OPTARG
        else
	           echo ":( Wrong syntax for minimal input contig/scaffold size. Exiting!" >&2
  	  	fi
    ;;
	*) echo "$USAGE" >&2
		exit 1
	;;
esac
done
shift $(( OPTIND-1 ))

## HANDLE ARGUMENTS

[ -z $1 ] || [ -z $2 ] && echo >&2 ":( Some input seems missing." && echo >&2 "$USAGE" && exit 1

current_cprops=$1
current_asm=$2

[ -z ${SIZE} ] && echo >&2 ":| Warning: no size limit was listed. Will put back all singletons without preferential alternative location." && SIZE=`awk '\$3>max{max=\$3}END{print max+1}' ${current_cprops}`

## HANDLE DEPENDENCIES

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

## MAIN FUNCTION

id=`basename ${current_cprops} .cprops`

# dump new asm and cprops
	
	awk -v size=${SIZE} -f ${pipeline}/seal/build-sealed-asm.awk ${current_cprops} ${current_asm} 2>${id}.sealed.cprops 1>${id}.sealed.asm
