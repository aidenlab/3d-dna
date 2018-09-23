#!/bin/bash

## Helper script to generate rainbow tracks illustrating liftover
## Written by: OD olga.dudchenko@bcm.edu

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

bin=500000;

USAGE="
*****************************
	./make-rainbow-tracks.sh -b <bin> -c <ref_chrom_number> ref.assembly target.assembly
*****************************
"

## HANDLE OPTIONS

while getopts "hc:b:" opt; do
case $opt in
    h) echo "$USAGE" >&1
        exit 0
    ;;
    c) chrom=$OPTARG 
    ;;
	b) bin=$OPTARG
    ;;
    *) echo "$USAGE" >&2
        exit 1
    ;;
esac
done

shift $(( OPTIND-1 ))

## HANDLE ARGUMENTS

[ -z $1 ] || [ -z $2 ] && echo >&2 ":( Some input seems missing." && echo >&2 "$USAGE" && exit 1

awk -f ${pipeline}/utils/convert-assembly-to-cprops-and-asm.awk $1
ref_cprops=${1::-9}.cprops
ref_asm=${1::-9}.asm
awk -f ${pipeline}/utils/convert-assembly-to-cprops-and-asm.awk $2
target_cprops=${2::-9}.cprops
target_asm=${2::-9}.asm

[ -z $chrom ] && chrom=`wc -l<${ref_asm}` && echo >&2 ":| Warning: Number of chromosomes in the reference genome is not listed. Coloring all $chrom scaffolds"

awk -v bin="$bin" -v chrom="$chrom" 'FILENAME==ARGV[1]{len[$2]=$3;next}FNR<=chrom{gsub("-","");for(i=1;i<=NF;i++){c+=len[$i]}}END{OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2", "color","id","X1","X2","Y1","Y2"; s=0; while((s+bin)<=c){print "assembly", s, s+bin, "assembly", s, s+bin, "0,0,0", s/c, s, s+bin, s, s+bin; s+=bin}; print "assembly", s, c, "assembly", s, c, "0,0,0", 1, s, c, s, c}' ${ref_cprops} ${ref_asm} | awk -f ${pipeline}/supp/lift-asm-annotations-to-input-annotations-w-split.awk ${ref_cprops} <(awk 'FILENAME==ARGV[1]{include[$2]=1;next}{print; gsub("-",""); for(i=1;i<=NF;i++){delete include[$i]}}END{for(i in include){print i}}' ${ref_cprops} ${ref_asm}) - | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${ref_cprops} <(awk '{print $2}' ${ref_cprops}) - | awk -f ${pipeline}/lift/lift-asm-annotations-to-input-annotations.awk ${target_cprops} <(awk '{print $2}' ${target_cprops}) - | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${target_cprops} ${target_asm} -

rm ${ref_cprops} ${ref_asm} ${target_cprops} ${target_asm}
