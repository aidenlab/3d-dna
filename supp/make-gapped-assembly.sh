#!/bin/bash

#### Description: Generates a gapped assembly file from .rawchrom.assembly and (original) fasta.
#### Input: assembly, fasta
#### Output: _HiC.assembly.
#### Dependencies: make-gap-bed.awk, edit-cprops, edit-fasta, edit-asm scripts
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 11/13/2018

USAGE="
*****************************************************
	./make-gapped-assembly.sh <modified.assembly> <original-input-fasta-file>
*****************************************************
"

gap_size=500

## HANDLE OPTIONS
while getopts "g:h" opt; do
case $opt in
	h) echo "$USAGE"
		exit 0
	;;
	g) 	re='^[0-9]+$'
		if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
		  echo "... -g flag was triggered, making gap size between scaffolded draft sequences to be equal to $OPTARG." >&1
		  gap_size=$OPTARG
		else
		  echo ":( Wrong syntax for default gap size parameter value. Using default gap_size=${gap_size}!" >&2
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

assembly=$1 # rawchrom
fasta=$2 # original

if [ ! -f ${assembly} ] || [ ! -f ${fasta} ]; then
	echo "Files corresponding to original input not found. Exiting!" >&2
	exit 1
fi

## HANDLE DEPENDENCIES

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

## remove suffixes

prefix=`basename ${assembly} .assembly`
#prefix=`basename ${prefix} .review`
# prefix=`basename ${prefix} .final`
# prefix=`basename ${prefix} .rawchrom`

echo $prefix


## reconstruct original assembly (do I need it?)
#${pipeline}/utils/reconstruct-assembly.sh ${assembly}


awk -f ${pipeline}/utils/convert-assembly-to-cprops-and-asm.awk ${assembly}

cprops=${prefix}.cprops
asm=${prefix}.asm

awk -f ${pipeline}/edit/edit-fasta-according-to-new-cprops.awk ${cprops} ${fasta} | awk -f ${pipeline}/utils/make-gap-bed.awk | awk 'BEGIN{OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2","color", "id", "X1", "X2", "Y1", "Y2"}FILENAME==ARGV[1]{len[$1]=$3;next}$2==0{print $1, $2, $3, $1, $2, $3, "0,0,0", "gap", $2, $3, $2, $3; next}$3==len[$1]{print $1, $2, $3, $1, $2, $3, "0,0,0", "gap", $2, $3, $2, $3}' ${cprops} - | awk -v label1=":::overhang_" -v label2=":::gap" -f ${pipeline}/edit/edit-cprops-according-to-annotations.awk - ${cprops} > ${prefix}.no_overhangs.cprops

bash ${pipeline}/edit/edit-asm-according-to-new-cprops.sh ${prefix}.no_overhangs.cprops ${cprops} ${asm} | awk 'FILENAME==ARGV[1]{if($1~/:::gap/){skip[$2]=1; skip[-$2]=1}; next}{str=""; for(i=1;i<=NF;i++){if(! skip[$i]){str=str" "$i}}; if(str!=""){print substr(str,2)}}' ${prefix}.no_overhangs.cprops - > ${prefix}.no_overhangs.asm
		
		cprops=${prefix}.no_overhangs.cprops
		asm=${prefix}.no_overhangs.asm
		
		# riffle
		echo "...adding gaps"
		
		gap_id=`awk 'END{print \$2+1}' ${cprops}`
	
		awk -v riffle=${gap_id} -f ${pipeline}/finalize/riffle-asm.awk ${asm} > temp.asm
	
		cp ${cprops} temp.cprops
		echo "hic_gap_${gap_id} ${gap_id} ${gap_size}" >> temp.cprops
		
		awk -v disable_checks=1 -f ${pipeline}/utils/convert-cprops-and-asm-to-assembly.awk temp.cprops temp.asm

		rm temp.cprops temp.asm ${cprops} ${asm} ${prefix}.cprops ${prefix}.asm		
		
		prefix=`basename ${prefix} .review`
		prefix=`basename ${prefix} .final`
		prefix=`basename ${prefix} .rawchrom`
		
		mv temp.assembly ${prefix}"_HiC.assembly"
