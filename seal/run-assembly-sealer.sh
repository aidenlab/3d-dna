#!/bin/bash
## Wrapper script to analyze the assembly and place fragments from among the small scaffolds back into the assembly
## NOTE: now is run after splitting but for diploid pipeline could be run for for tiled assembly
## NOTE: Relies on standard annotations :::fragment_ and :::debris!
## NOTE: Probably should be done after chrom splitting to avoid a slim chance that a misassembled contig/scaffold spanned chromosomes, and the chromosomes ended up in the order and orientation to fit the edges of the contig/scaffold

USAGE="
**********************************
	
	./run-assembly-sealer.sh -i input_size -b bundle_size <current_assembly>
	
**********************************
"

## HANDLE OPTIONS

while getopts "i:b:h" opt; do
case $opt in
	h) echo "$USAGE"
		exit 0
	;;
	i) 	re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -ge 0 ]]; then
	          echo "... -i flag was triggered, will attempt to place back only debris contigs/scaffolds and those less than $OPTARG" >&1
	          SIZE=$OPTARG
        else
	           echo ":( Wrong syntax for minimal input contig/scaffold size to consider as candidates for false positive removal. Exiting!" >&2
  	  	fi
    ;;
    b)	re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -ge 0 ]]; then
	          echo "... -b flag was triggered, will attempt to group all contigs/scaffolds smaller than $OPTARG at the end of the assembly" >&1
	          BUNDLE=$OPTARG
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

[ -z $1 ] && echo >&2 ":( Some input seems missing." && echo >&2 "$USAGE" && exit 1

current_assembly=$1

# if no explicit requirement re bundle size is made use sensitivity size as this is most typical scenario
([ -z ${BUNDLE} ] && [ ! -z ${SIZE} ]) && echo >&2 ":| Warning: no explicit bundle size was listed. Will use the same one as listed for false positive size threshold: this is the most typical scenario." && BUNDLE=${SIZE}
# if no input size is passed try to place back everything
[ -z ${SIZE} ] && echo >&2 ":| Warning: no size limit was listed. Will put back all singletons without preferential alternative location, regardless of size." && SIZE=`awk '\$0!~/^>/{next}\$3>max{max=\$3}END{print max+1}' ${current_assembly}`
# if no size and no bundle is passed do not bundle
[ -z ${BUNDLE} ] && echo >&2 ":| Warning: no size limit and no bundle size was listed. Will not bundle." && BUNDLE=0 

## HANDLE DEPENDENCIES

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

## MAIN FUNCTION

id=`basename ${current_assembly} .assembly`

# dump new assembly
	
awk -v size=${SIZE} -v bundle=${BUNDLE} -f ${pipeline}/seal/seal-assembly.awk ${current_assembly} >${id}.sealed.assembly
