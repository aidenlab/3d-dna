#!/bin/bash

## Part of the merge segment of the diploid pipeline.
## Requires LASTZ in path! Requires PARALLEL in path!
## TODO: Add scaling for annotations? Or fully switch to bpe? Add warning!!
## Written by: Olga Dudchenko

USAGE="
*****************************
	./run-pairwise-alignment.sh [options] <path_to_assembly> <path_to_fasta>
*****************************
"
# Defaults:
merger_search_band=3000000
merger_alignment_score=50000000
merger_alignment_identity=20
merger_alignment_length=20000
merger_lastz_options=\"--gfextend\ --gapped\ --ydrop=300000\ --hspthresh=50000\ ‑‑allocate:traceback=2.0G\ --noytrim\ --gappedthresh=1000000\ --step=5000\ --chain\"


#echo ${merger_lastz_options}

#--chain=200,200\ 
#--strand=0
#--masking=1

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

[ -z $1 ] || [ -z $2 ] && echo >&2 ":( Some input seems missing." && echo >&2 "$USAGE" && exit 1

assembly=$1
fasta=$2

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

awk -v band=${merger_search_band} '$0~/^>/{len[$2]=$3; oname[$2]=substr($1,2); next}{gsub("-",""); for(i=1; i<NF; i++){shift=0; str=""; k=i+1; while( shift<=band && k<=NF ){str=str","oname[$k]; shift+=len[$k]; k++}; if(str!=""){print oname[$i], substr(str,2)};}}' $assembly > joblist.txt

echo "chr1	x1	x2	chr2	y1	y2	color	id	X1	X2	Y1	Y2" > overlaps_2D_input_raw.txt

# could have used lastz indexing but hard to find tools..
[ -f fasta.index ] && rm fasta.index
#mkfifo fasta.index
awk -f ${pipeline}/utils/index-lines-in-fasta.awk ${fasta} > fasta.index #&

if [ $parallel == "true" ]; then
	parallel --colsep=" " --will-cite --jobs 80% -a joblist.txt "lastz <(bash ${pipeline}/finalize/print-fasta-subset.sh -i fasta.index {1} ${fasta})[unmask] <(bash ${pipeline}/finalize/print-fasta-subset.sh -i fasta.index {2} ${fasta})[unmask] ${merger_lastz_options}" | awk -v skip_header=1 -f ${pipeline}/merge/extract-all-stanzas-and-convert-to-input-annotations.awk - >> overlaps_2D_input_raw.txt 
else
	while read target query
	do
		lastz <(bash ${pipeline}/finalize/print-fasta-subset.sh -i fasta.index "$target" ${fasta})[unmask] <(bash ${pipeline}/finalize/print-fasta-subset.sh -i fasta.index "$query" ${fasta})[unmask] ${merger_lastz_options} | awk -v skip_header=1 -f ${pipeline}/merge/extract-all-stanzas-and-convert-to-input-annotations.awk - >> overlaps_2D_input_raw.txt
	done < joblist.txt
fi

rm fasta.index

## we expect no more than one overlap per pair of input sequences: nearby and overlapping pairs should be merged

echo "chr1	x1	x2	chr2	y1	y2	color	id	X1	X2	Y1	Y2" > overlaps_2D_input.txt

awk 'NR>1{print $0"\t"substr($8,2)}' overlaps_2D_input_raw.txt | sort -k 1,1 -k 4,4 -k 13,13nr | awk 'BEGIN{FS="\t"; OFS="\t"}($1!=target||$4!=query){if(target){print target, target_start, target_end, query, query_start, query_end, "0,0,0", score, target_start, target_end, query_start, query_end}; target=$1; query=$4; target_start=$2; target_end=$3; query_start=$5; query_end=$6; score=$8; next}($8>0&&target_end<$3&&query_end<$6){target_end=$3;query_end=$6}($8>0&&target_start>$2&&query_start>$5){target_start=$2;query_start=$5}($8<0&&target_end<$3&&query_start>$5){target_end=$3;query_start=$5}($8<0&&target_start>$2&&query_end<$6){target_start=$2;query_end=$6}END{print target, target_start, target_end, query, query_start, query_end, "0,0,0", score, target_start, target_end, query_start, query_end}' >> overlaps_2D_input.txt

## temp

#awk -v min_overlap=${merger_alignment_length} -v identity=${merger_alignment_identity} -v max_k=${merger_alignment_score} '{if($6-$5<$9-$8){overlap=($6-$5)}else{overlap=($9-$8)}; tol=.5*overlap}(overlap>min_overlap)&&($1>identity*overlap||$1>max_k)&&(($10==0&&($5<tol||$8<tol)&&($2-$6<tol||$3-$9<tol))||($10==1&&($5<tol||$3-$9<tol)&&($2-$6<tol||$8<tol))){print}' alignments.txt | awk 'BEGIN{OFS="\t"}{for(i=1; i<=NF-3; i++){$i=$(i+3)}; $2--; $5--; $8=$7; $7="0,0,0"; $9=$2; $10=$3; $11=$5; $12=$6}1' >> overlaps_2D_input.txt

awk -f ${pipeline}/lift/lift-input-annotations-to-assembly-annotations.awk $assembly overlaps_2D_input.txt > overlaps_2D_asm.txt


