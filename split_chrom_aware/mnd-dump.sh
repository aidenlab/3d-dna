#!/bin/bash

#### Description: Script to dump matrix data from mnd.
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu.

#TODO: Add usage

## Set defaults

mapq=1

## Handle options

while getopts "q:" opt; do
case $opt in
    q)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            mapq=$OPTARG
        else
            echo ":( Wrong syntax for mapping quality. Using the default value ${mapq}" >&2
        fi
    ;;
    *) echo "$USAGE" >&2
        exit 1
    ;;
esac
done

shift $(( OPTIND-1 ))

## Handle arguments
mnd=$1
bin_size=$2

[ -s ${mnd} ] && [ ! -z ${bin_size} ] || { echo >&2 ":( Not sure how to parse your input or input files not found at intended locations. Exiting!" && exit 1 ; }

## Main function

##	GNU Parallel Dependency
parallel="false"
if hash parallel 2>/dev/null; then
        ver=`parallel --version | awk 'NR==1{print \$3}'`
        [ $ver -ge 20150322 ] && parallel="true"
fi
[ $parallel == "false" ] && echo >&2 ":| WARNING: GNU Parallel version 20150322 or later not installed. We highly recommend to install it to increase performance. Starting pipeline without parallelization!"

cmd="awk -v b=${bin_size} -v q=${mapq} '(\$9>=q)&&(\$12>=q)&&(\$4!=\$8){if(\$3<=\$7){c[b*int(\$3/b)\" \"b*int(\$7/b)]+=1}else{c[b*int(\$7/b)\" \"b*int(\$3/b)]+=1}}END{for(i in c){print i, c[i]}}'"

[ $parallel == "true" ] && parallel --will-cite -a ${mnd} --pipepart -j 80% --block 1G ${cmd} | awk -v b=${bin_size} -v OFS='\t' '{c[$1"\t"$2]+=$3}END{for(i in c){print i, c[i]}}'

[ $parallel == "false" ] && { eval ${cmd}" "${mnd} | awk '{print $1"\t"$2"\t"$3}' -; }
