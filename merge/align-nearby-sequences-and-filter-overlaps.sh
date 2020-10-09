#!/bin/bash

## Part of the merge segment of the diploid pipeline.
## Requires LASTZ in path! Requires PARALLEL in path!
## TODO: Add scaling for annotations? Or fully switch to bpe? Add warning!!
## Written by: Olga Dudchenko

USAGE="
*****************************
	./run-pairwise-alignment.sh [options] <path_to_cprops> <path_to_asm> <path_to_faSplitDir>
*****************************
"
# Defaults:
merger_search_band=3000000
merger_alignment_score=50000000
merger_alignment_identity=20
merger_alignment_length=20000
merger_lastz_options=\"--gfextend\ --gapped\ --chain=200,200\"

## HANDLE OPTIONS

while getopts "b:s:i:l:o:h" opt; do
case $opt in
	b)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
            echo ":) -b flag was triggered, assuming ${OPTARG} bp-wide band for alternative haplotype detection." >&1
            merger_search_band=$OPTARG
        else
            echo ":( Wrong syntax for band size. Using the default value ${band_size}." >&2
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

# HANDLE DEPENDENCIES
pipeline=`cd "$( dirname $0)" && cd .. && pwd`

## LASTZ dependency
lastz="false"
if hash lastz 2>/dev/null; then
	lastz="true"
fi
[[ $lastz == "false" ]] && echo >&2 ":( LASTZ not installed or not in the path. The merge section of the pipeline will not run. Exiting!" && exit 1

##	GNU Parallel Dependency
parallel="false"
if hash parallel 2>/dev/null; then
        ver=`parallel --version | awk 'NR==1{print \$3}'`
        [ $ver -ge 20150322 ] && parallel="true"
fi
[ $parallel == "false" ] && echo ":| WARNING: GNU Parallel version 20150322 or later not installed. We highly recommend to install it to increase performance. Starting pipeline without parallelization!" >&2

# remove comments from LASTZ options
merger_lastz_options="${merger_lastz_options%\"}"
merger_lastz_options="${merger_lastz_options#\"}"

# read in the lengths as associative array

declare -A len
while read -r -a array
do 
  	len["${array[1]}"]="${array[2]}"
done < $cprops

# cleanup

[ -f alignments.txt ] &&  rm alignments.txt
[ -f joblist.txt ] &&  rm joblist.txt

# read asm line by line

while read -r line
do
	arr=($line)	
	for ((index=0; index <= ${#arr[@]}; index++)); do
		shift=0
		k=$((index+1))
		
		while [ $shift -le ${merger_search_band} ] && [ $k -lt ${#arr[@]} ]; do
			## create a job to align $index to $k
			echo ${arr[index]}" "${arr[k]} >> joblist.txt
			## count band from lower end of the contig/scaffold
			shift=$((shift + ${len[${arr[k]}]}))
			k=$((k+1))
		done	
	done

done < <(awk '{gsub("-",""); print}' $asm)

if [ $parallel == "true" ]; then
	parallel --colsep=" " --will-cite --jobs 80% -a joblist.txt "lastz faSplit/{1}.fa faSplit/{2}.fa ${merger_lastz_options} | awk -f ${pipeline}/merge/extract-total-stanza.awk" > alignments.txt
else
	while read oname1 oname2
	do
		lastz faSplit/${oname1}.fa faSplit/${oname2}.fa ${merger_lastz_options} | awk -f ${pipeline}/merge/extract-total-stanza.awk >> alignments.txt
	done < joblist.txt
fi 

echo "chr1	x1	x2	chr2	y1	y2	color	id	X1	X2	Y1	Y2" > overlaps_2D_input.txt

awk -v min_overlap=${merger_alignment_length} -v identity=${merger_alignment_identity} -v max_k=${merger_alignment_score} '{if($6-$5<$9-$8){overlap=($6-$5)}else{overlap=($9-$8)}; tol=.5*overlap}(overlap>min_overlap)&&($1>identity*overlap||$1>max_k)&&(($10==0&&($5<tol||$8<tol)&&($2-$6<tol||$3-$9<tol))||($10==1&&($5<tol||$3-$9<tol)&&($2-$6<tol||$8<tol))){print}' alignments.txt | awk 'BEGIN{OFS="\t"}{for(i=1; i<=NF-3; i++){$i=$(i+3)}; $2--; $5--; $8=$7; $7="0,0,0"; $9=$2; $10=$3; $11=$5; $12=$6}1' >> overlaps_2D_input.txt

awk -f ${pipeline}/supp/lift-input-annotations-to-asm-annotations.awk ${cprops} ${asm} overlaps_2D_input.txt > overlaps_2D_asm.txt
