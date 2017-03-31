#!/bin/bash

## Part of the merge segment of the diploid pipeline.
## Requires LASTZ in path!
## TODO: Add scaling for annotations. Add warning!!
## Written by: Olga Dudchenko

USAGE="
*****************************
	./run-pairwise-alignment.sh <path_to_cprops> <path_to_asm> <path_to_faSplitDir>
*****************************
"
# Defaults:
band_size=1000000	# band in which to search for alternative haplotypes

min_overlap=10000	# criterion
identity=25
max_k=2500000 

lastz_opt="--gfextend --gapped --chain"	# lastz opt string


## HANDLE OPTIONS

while getopts "b:h" opt; do
case $opt in
	b)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
            echo ":) -b flag was triggered, assuming $OPTARG wide band alternative haplotype detection" >&1
            band_size=$OPTARG
        else
            echo ":( Wrong syntax for band size. Using the default value ${band_size}" >&2
        fi
    ;;
    h) echo "$USAGE" >&1
        exit 0
    ;;
    *) echo "$USAGE" >&2
        exit 1
    ;;
esac
done

shift $(( OPTIND-1 ))


## HANDLE ARGUMENTS

[ -z $1 ] || [ -z $2 ] || [ -z $3 ] && echo >&2 ":( Some input seems missing." && echo >&2 "$USAGE" && exit 1

cprops=$1
asm=$2
faSplit=$3

# Handle dependencies
pipeline=`cd "$( dirname $0)" && cd .. && pwd`

##	GNU Parallel Dependency
parallel="false"
if hash parallel 2>/dev/null; then
        ver=`parallel --version | awk 'NR==1{print \$3}'`
        [ $ver -ge 20150322 ] && parallel="true"
fi

[ $parallel == "false" ] && echo >&2 ":| WARNING: GNU Parallel version 20150322 or later not installed. We highly recommend to install it to increase performance. Starting pipeline without parallelization!"


# read in the lengths as associative array

declare -A len
while read -r -a array
do 
  	len["${array[1]}"]="${array[2]}"
done < $cprops

# read asm line by line

#echo "chr1	x1	x2	chr2	y1	y2	color	id	X1	X2	Y1	Y2" > overlaps_2D_input.txt

[ -f alignments.txt ] &&  rm alignments.txt
[ $parallel == "true" ] && [ -f joblist.txt ] &&  rm joblist.txt

while read -r line
do
	arr=($line)
	
	for ((index=0; index <= ${#arr[@]}; index++)); do

		shift=0
		k=$((index+1))
		
		while [ $shift -le $band_size ] && [ $k -lt ${#arr[@]} ]; do
			## align $index to $k
 			[ $parallel == "true" ] && echo ${arr[index]}" "${arr[k]} >> joblist.txt
 			[ $parallel == "false" ] && { align=$(lastz faSplit/${arr[index]}.fa faSplit/${arr[k]}.fa ${lastz_opt} | awk -f ${pipeline}/merge/extract-total-stanza.awk); [[ $align == "" ]] || echo $align >> alignments.txt; }

			## count band from lower end of the contig/scaffold
			shift=$((shift + ${len[${arr[k]}]}))
			k=$((k+1))
		done
	
	done

done < <(awk '{gsub("-",""); print}' $asm)

[ $parallel == "true" ] && { parallel --colsep=" " --will-cite --jobs 80% -a joblist.txt "lastz faSplit/{1}.fa faSplit/{2}.fa ${lastz_opt} | awk -f ${pipeline}/merge/extract-total-stanza.awk" > alignments.txt; }

echo "chr1	x1	x2	chr2	y1	y2	color	id	X1	X2	Y1	Y2" > overlaps_2D_input.txt

awk -v min_overlap=${min_overlap} -v identity=${identity} -v max_k=${max_k} '{if($6-$5<$9-$8){overlap=($6-$5)}else{overlap=($9-$8)}; tol=.5*overlap}(overlap>min_overlap)&&($1>identity*overlap||$1>max_k)&&(($10==0&&($5<tol||$8<tol)&&($2-$6<tol||$3-$9<tol))||($10==1&&($5<tol||$3-$9<tol)&&($2-$6<tol||$8<tol))){print}' alignments.txt | awk 'BEGIN{OFS="\t"}{for(i=1; i<=NF-3; i++){$i=$(i+3)}; $2--; $5--; $8=$7; $7="0,0,0"; $9=$2; $10=$3; $11=$5; $12=$6}1' >> overlaps_2D_input.txt

# rectangular
awk -f ${pipeline}/supp/lift-input-annotations-to-asm-annotations.awk ${cprops} ${asm} overlaps_2D_input.txt > overlaps_2D_asm.txt
