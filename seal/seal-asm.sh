#!/bin/bash
## Wrapper script to analyze the assembly and place fragments from among the small scaffolds back into the assembly
## NOTE: now is run after splitting but for diploid pipeline could be run for for tiled assembly 
## NOTE: Relies on standard annotations :::fragment_ and :::debris!
## NOTE: Probably should be done after chrom splitting to avoid a slim chance that a misassembled contig/scaffold spanned chromosomes, and the chromosomes ended up in the order and orientation to fit the edges of the contig/scaffold

USAGE="
**********************************
	
	./seal-asm.sh <orig_cprops> <current_cprops> <current_asm>
	
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

[ -z $1 ] || [ -z $2 ] || [ -z $3 ] && echo >&2 ":( Some input seems missing." && echo >&2 "$USAGE" && exit 1

orig_cprops=$1
current_cprops=$2
current_asm=$3

[ -z ${SIZE} ] && echo >&2 ":| Warning: no size limit was listed. Will put back all singletons without preferential alternative location." && SIZE=`awk '\$3>max{max=\$3}END{print max+1}' ${current_cprops}`

## HANDLE DEPENDENCIES

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

## MAIN FUNCTION

id=`basename ${current_cprops} .cprops`

#	1) Dump new asm. Create lits of annotations to be deleted (filename1) and to be merged (filename2) in the process

	touch "h."${id}".to_delete.txt"
	touch "h."${id}".to_merge.txt"
	
	awk -v filename1="h."${id}".to_delete.txt" -v filename2="h."${id}".to_merge.txt" -v size=${SIZE} -f ${pipeline}/seal/build-sealed-asm.awk ${current_cprops} ${current_asm} > ${id}.sealed.asm

#	2) Convert annotations into the orig input coordinate system
		
	awk 'BEGIN{OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "id", "X1", "X2", "Y1", "Y2"}FILENAME==ARGV[1]{oname[$2]=$1; len[$2]=$3; next}{for(i=1; i<=NF; i++){if(oname[$i]~/:::debris/){print oname[$i], 0, len[$i], oname[$i], 0, len[$i], "0,0,0", "debris", 0, len[$i], 0, len[$i]}}}' ${current_cprops} "h."${id}".to_delete.txt" | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${current_cprops} <(awk '{print $2}' ${current_cprops}) - | awk -f ${pipeline}/lift/lift-asm-annotations-to-input-annotations.awk ${orig_cprops} <(awk '{print $2}' ${orig_cprops}) - > "h."${id}".to_delete.txt.tmp" 
	
	awk 'BEGIN{OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "id", "X1", "X2", "Y1", "Y2"}FILENAME==ARGV[1]{oname[$2]=$1; len[$2]=$3; next}{for(i=1; i<=NF; i++){if(oname[$i]~/:::debris/){print oname[$i], 0, len[$i], oname[$i], 0, len[$i], "0,0,0", "debris", 0, len[$i], 0, len[$i]}}}' ${current_cprops} "h."${id}".to_merge.txt" | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${current_cprops} <(awk '{print $2}' ${current_cprops}) - | awk -f ${pipeline}/lift/lift-asm-annotations-to-input-annotations.awk ${orig_cprops} <(awk '{print $2}' ${orig_cprops}) - | awk 'NR>1' >> "h."${id}".to_delete.txt.tmp"
	
	awk 'BEGIN{OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "id", "X1", "X2", "Y1", "Y2"}FILENAME==ARGV[1]{oname[$2]=$1; len[$2]=$3; next}{print oname[$1], 0, 1, oname[$1], 0, 1, "0,0,0", "debris", 0,1,0,1; print oname[$NF], len[$NF]-1, len[$NF], oname[$NF], len[$NF]-1, len[$NF], "0,0,0", "debris", len[$NF]-1,len[$NF],len[$NF]-1,len[$NF]}' ${current_cprops} "h."${id}".to_merge.txt" | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${current_cprops} <(awk '{print $2}' ${current_cprops}) - | awk -f ${pipeline}/lift/lift-asm-annotations-to-input-annotations.awk ${orig_cprops} <(awk '{print $2}' ${orig_cprops}) - | awk 'BEGIN{OFS="\t"}NR==1{print; next}NR/2==int(NR/2){start=$2; next}{$2=start; $5=start; $9=start; $11=start; print}' > "h."${id}".to_merge.txt.tmp"

#	2) Reconstruct annotations from current_cprops and make the requested changes
	
	awk 'BEGIN{OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "id", "X1", "X2", "Y1", "Y2"}$1~/:::debris/{print $1, 0, $3, $1, 0, $3, "0,0,0", "debris", 0, $3, 0, $3}' ${current_cprops} | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${current_cprops} <(awk '{print $2}' ${current_cprops}) - | awk -f ${pipeline}/lift/lift-asm-annotations-to-input-annotations.awk ${orig_cprops} <(awk '{print $2}' ${orig_cprops}) - | awk 'FILENAME==ARGV[1]{skip[$0]=1; next}(!skip[$0])&&(FNR>1){print}' "h."${id}".to_delete.txt.tmp" - > ${id}_edit_annotations_sealed.txt
	
	awk 'NR>1' "h."${id}".to_merge.txt.tmp" ${id}_edit_annotations_sealed.txt | sort -k 1,1 -k 2,2n -k 3,3n > ${id}_edit_annotations_sealed.txt.tmp
	
	#add header
	awk 'BEGIN{OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "id", "X1", "X2", "Y1", "Y2"}' > ${id}_edit_annotations_sealed.txt
	
	cat ${id}_edit_annotations_sealed.txt.tmp >> ${id}_edit_annotations_sealed.txt

#	3) Edit cprops according to annotations
	awk -f ${pipeline}/edit/edit-cprops-according-to-annotations.awk ${id}_edit_annotations_sealed.txt ${orig_cprops} > ${id}.sealed.cprops
	
#	4) Clean up
	rm ${id}_edit_annotations_sealed.txt.tmp
	rm "h."${id}".to_merge.txt.tmp" "h."${id}".to_merge.txt"
	rm "h."${id}".to_delete.txt.tmp" "h."${id}".to_delete.txt"
	
	
