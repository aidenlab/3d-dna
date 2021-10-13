#!/bin/bash

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

USAGE="

	./edit-and-lift-input-annotations-to-assembly-annotations.sh -b <bundle-size> <path_to_assembly> <path_to_annotations>

"

# handle options

while getopts "hb:" opt; do
    case "$opt" in
    h)
        echo "$USAGE"
        exit 0
        ;;
    b)  # add safeguards
    	bundle_size=$OPTARG	
        ;;
    esac
done

shift $((OPTIND-1))


# handle arguments

assembly=$1
annotations=$2

# 1) reconstruct original assembly

bundled_assembly=${assembly}

if [ ! -z ${bundle_size} ]; then
	awk -v input_size=${bundle_size} -f ${pipeline}/utils/bundle-unattempted.awk ${assembly}
	bundled_assembly=`basename ${assembly} .assembly`".bundled.assembly"
fi

# 2) lift input annotations to original assembly, lift original to edit, lift edit to edit rearranged

 awk -f ${pipeline}/lift/lift-input-annotations-to-assembly-annotations.awk <(awk '$0~/>/{split($1,a,":::fragment_||:::debris"); if(!(a[1] in id)){counter++; id[a[1]]=1; name[counter]=a[1]}; len[a[1]]+=$3}END{for(i=1;i<=length(name);i++){print name[i], i, len[name[i]]}; for(i=1;i<=length(name); i++){print i}}' ${assembly}) ${annotations} | awk -f ${pipeline}/lift/lift-assembly-annotations-to-input-annotations.awk <(awk '$0~/^>/{counter++; print}END{for(i=1;i<=counter;i++){print i}}' ${bundled_assembly}) - | awk -f ${pipeline}/lift/lift-input-annotations-to-assembly-annotations.awk ${bundled_assembly} -

if [ ! -z ${bundle_size} ]; then
	rm `basename ${assembly} .assembly`".bundled.assembly"
	rm `basename ${assembly} .assembly`".unattempted.assembly"
fi