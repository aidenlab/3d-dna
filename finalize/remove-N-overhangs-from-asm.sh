#!/bin/bash

#### Description: Script to create an equivalent to a given assembly output (cprops asm fasta) but without N overhangs in the input contigs/scaffolds.
#### Usage: bash remove-N-overhangs-from-asm.sh <(modified) input cprops> <(modified) input asm> <(modified) input-fasta-file>
#### Input: cprops, asm, fasta
#### Output: "no_overhangs" cprops, asm, fasta.
#### Dependencies: make-gap-bed.awk, edit-cprops, edit-fasta, edit-asm scripts
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 01/30/2017

USAGE="
*****************************************************
	./remove-N-overhangs-from-asm.sh <(edited) input cprops> <(edited) input asm> <(edited) input-fasta-file>
*****************************************************
"
## HANDLE OPTIONS
while getopts "h" opt; do
case $opt in
	h) echo "$USAGE"
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

orig_cprops=$1
orig_asm=$2
orig_fasta=$3

if [ ! -f ${orig_cprops} ] || [ ! -f ${orig_asm} ] || [ ! -f ${orig_fasta} ]; then
	echo "Files corresponding to original input not found. Exiting!" >&2
	exit 1
fi

## HANDLE DEPENDENCIES

pipeline=`cd "$( dirname $0)" && cd .. && pwd`
make_gap_bed=${pipeline}/utils/make-gap-bed.awk
edit_cprops=${pipeline}/edit/edit-cprops-according-to-annotations.awk
edit_asm=${pipeline}/edit/edit-asm-according-to-new-cprops.sh
#edit_fasta=${pipeline}/edit/edit-fasta-according-to-new-cprops.awk
edit_fasta=${pipeline}/edit/edit-fasta-according-to-new-cprops.py

prefix=`basename ${orig_cprops} .cprops`

awk -f ${make_gap_bed} ${orig_fasta} | \
	awk 'BEGIN{OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2","color", "id", "X1", "X2", "Y1", "Y2"}FILENAME==ARGV[1]{len[$1]=$3;next}$2==0{print $1, $2, $3, $1, $2, $3, "0,0,0", "gap", $2, $3, $2, $3; next}$3==len[$1]{print $1, $2, $3, $1, $2, $3, "0,0,0", "gap", $2, $3, $2, $3}' ${orig_cprops} - | \
	awk -v label1=":::overhang_" -v label2=":::gap" -f ${edit_cprops} - ${orig_cprops} > ${prefix}.no_overhangs.cprops

new_cprops=${prefix}.no_overhangs.cprops

prefix=`basename ${orig_asm} .asm`

bash ${edit_asm} ${new_cprops} ${orig_cprops} ${orig_asm} | \
	awk 'FILENAME==ARGV[1]{if($1~/:::gap/){skip[$2]=1; skip[-$2]=1}; next}{str=""; for(i=1;i<=NF;i++){if(! skip[$i]){str=str" "$i}}; if(str!=""){print substr(str,2)}}' ${new_cprops} - > ${prefix}.no_overhangs.asm

prefix=`basename ${orig_fasta} .fa`
prefix=`basename ${prefix} .fna`
prefix=`basename ${prefix} .fasta`

#awk -v label1=":::overhang_" -v label2=":::gap" -f ${edit_fasta} ${new_cprops} ${orig_fasta} | awk '$0~/>/{test=1; if($0~/:::gap/){test=0}}test{print}' > ${prefix}.no_overhangs.fasta
python3  ${edit_fasta} --label overhang ${new_cprops} ${orig_fasta} | seqkit seq -w 80 | awk '$0~/>/{test=1; if($0~/:::gap/){test=0}}test{print}' > ${prefix}.no_overhangs.fasta

