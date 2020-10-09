#!/bin/bash
# Script to make an equivalent of the old asm file in terms of new cprops if contigs/scaffolds were edited
# Written by OD

new_cprops=$1
old_cprops=$2
old_asm=$3

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

## TODO: SAFEGUARDS

((awk 'BEGIN{OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "id", "X1", "X2", "Y1", "Y2"}{print $1, 0, $3, $1, 0, $3, "0,0,0", $2, 0, $3, 0, $3}' ${old_cprops} | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${old_cprops} <(awk '{print $2}' ${old_cprops}) - | awk 'NR>1{print $2, $3, $8, 2}') && (awk 'BEGIN{OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "id", "X1", "X2", "Y1", "Y2"}{print $1, 0, $3, $1, 0, $3, "0,0,0", $2, 0, $3, 0, $3}' ${new_cprops} | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${new_cprops} <(awk '{print $2}' ${new_cprops}) - | awk 'NR>1{print $2, $3, $8, 1}')) | sort -k 2,2n -k 4,4n | awk '$4==2{print $3, substr(new_id,2); new_id=""; next}{new_id=new_id" "$3}' | awk 'function reverse(seq){n=split(seq,s); tmpstring = -s[n]; for(k=n-1; k>=1; k--){tmpstring = sprintf("%s %s",tmpstring, -s[k])};return tmpstring}FILENAME==ARGV[1]{str=$2; for(i=3;i<=NF;i++){str=str" "$i}; substitute[$1]=str; next}{str=""; for(i=1;i<=NF;i++){if($i<0){str=str" "reverse(substitute[-$i])}else{str=str" "substitute[$i]}}; print substr(str,2)}' - ${old_asm} 
