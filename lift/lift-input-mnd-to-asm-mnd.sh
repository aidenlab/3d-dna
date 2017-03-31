#!/bin/bash

## Handle options
scale=1
while getopts "s:" opt; do
case $opt in
    s)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            scale=$OPTARG
        else
            echo ":( Wrong syntax for scale quality. Using the default value ${scale}" >&2
        fi
    ;;
    *) echo "$USAGE" >&2
        exit 1
    ;;
esac
done

shift $(( OPTIND-1 ))

cprops=$1
asm=$2
mnd=$3

pipeline=`cd "$( dirname $0)" && cd .. && pwd`
lift_mnd_script=${pipeline}/lift/lift-input-mnd-to-asm-mnd.awk

##	GNU Parallel Dependency
parallel="false"
if hash parallel 2>/dev/null; then
        ver=`parallel --version | awk 'NR==1{print \$3}'`
        [ $ver -ge 20150322 ] && parallel="true"
fi

[ $parallel == "false" ] && echo >&2 ":| WARNING: GNU Parallel version 20150322 or later not installed. We highly recommend to install it to increase performance. Starting pipeline without parallelization!"

[ $parallel == "true" ] && parallel -a ${mnd} --pipepart --will-cite --jobs 80% --block 1G "awk -v scale=${scale} -f ${lift_mnd_script} ${cprops} ${asm} - "
[ $parallel == "false" ] && awk -v scale=${scale} -f ${lift_mnd_script} ${cprops} ${asm} ${mnd}
