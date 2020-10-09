#!/bin/bash

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

assembly=$1
annotations=$2


# 1) reconstruct original assembly

mkfifo ${assembly}"_original_assembly.tmp"
awk '$0~/>/{split($1,a,":::fragment_||:::debris"); if(!(a[1] in id)){counter++; id[a[1]]=1; name[counter]=a[1]}; len[a[1]]+=$3}END{for(i=1;i<=length(name);i++){print name[i], i, len[name[i]]}; for(i=1;i<=length(name); i++){print i}}' ${assembly} > ${assembly}"_original_assembly.tmp" &

# 2) lift input annotations to original assembly, lift original to edit, lift edit to edit rearranged

 awk -f ${pipeline}/lift/lift-input-annotations-to-assembly-annotations.awk ${assembly}"_original_assembly.tmp" $annotations | awk -f ${pipeline}/lift/lift-assembly-annotations-to-input-annotations.awk <(awk '$0~/^>/{counter++; print}END{for(i=1;i<=counter;i++){print i}}' ${assembly}) - | awk -f ${pipeline}/lift/lift-input-annotations-to-assembly-annotations.awk ${assembly} -
 
rm ${assembly}"_original_assembly.tmp"