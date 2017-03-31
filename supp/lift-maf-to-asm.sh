#!/bin/bash

## NOTE: Assumes that the genome that is being assembled is query in maf! TODO: do both ways

USAGE="
************************
./lift-maf-to-asm.sh -t <target_label_in_maf> -q <query_label_in_maf> <path_to_orig_cprops> <path_to_current_cprops> <path_to_current_asm> <path_to_orig_maf>
************************
"

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

## HANDLE OPTIONS

# Unprompted
bin=1000000

while getopts "t:q:rh" opt; do
case $opt in
h) echo "$USAGE"
exit 0
;;
t) id1=$OPTARG
	echo "...assuming target label in maf to be ${OPTARG}"
;;
q) id2=$OPTARG
	echo "...assuming query label in maf to be ${OPTARG}"
;;
r) rename=true
;;
*) echo "$USAGE" >&2
exit 1
;;
esac
done
shift $(( OPTIND-1 ))

[ -z $id1 ] || [ -z $id2 ] && echo >&2 "Please list maf labels -t and -q. Exiting!" && exit 1

original_cprops=$1
final_cprops=$2
final_asm=$3
maf_file=$4

[[ "$rename" == "true" ]] && rename_table=${pipeline}/data/genbank_to_supercontig.txt

prefix=`basename $maf_file .maf`
#echo "Running map for: "$prefix

echo "variableStep chrom=assembly span=""$bin" > "$prefix"".synteny_track.wig"

[ -z $rename_table ] && awk -v id1=${id1} -v id2=${id2} -f ${pipeline}/supp/maf-to-annotations.awk ${maf_file} | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${original_cprops} <(awk '{print $2}' ${original_cprops}) - | awk -f ${pipeline}/supp/lift-asm-annotations-to-input-annotations.awk ${final_cprops} <(awk '{print $2}' ${final_cprops}) - | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${final_cprops} ${final_asm} - | awk -v bin="$bin" '{count[int($2/bin)]+=1}END{for(i in count){print i*bin"\t"count[i]}}' - | sort -k 1,1n >> "$prefix"".synteny_track.wig"


[ -z $rename_table ] || awk -v id1=${id1} -v id2=${id2} -f ${pipeline}/supp/maf-to-annotations.awk ${maf_file} | awk 'BEGIN{OFS="\t"}FILENAME==ARGV[1]{genbankname[$2]=$1;next}{FS="\t"}{$1=genbankname[$1]; $4=genbankname[$4]}1' ${rename_table} - | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${original_cprops} <(awk '{print $2}' ${original_cprops}) - | awk -f ${pipeline}/supp/lift-asm-annotations-to-input-annotations.awk ${final_cprops} <(awk '{print $2}' ${final_cprops}) - | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${final_cprops} ${final_asm} - | awk -v bin="$bin" '{count[int($2/bin)]+=1}END{for(i in count){print i*bin"\t"count[i]}}' - | sort -k 1,1n >> "$prefix"".synteny_track.wig"



##
##
