#!/bin/bash
# Wrapper around gnu parallel for editing the mnd.
# ./edit-mnd-according-to-new-assembly.sh <new-assembly> <old-mnd>
# Written by: OD

##	GNU Parallel Dependency
parallel="false"
if hash parallel 2>/dev/null; then
        ver=`parallel --version | awk 'NR==1{print \$3}'`
        [ $ver -ge 20150322 ] && parallel="true"
fi
[ $parallel == "false" ] && echo >&2 ":| WARNING: GNU Parallel version 20150322 or later not installed. We highly recommend to install it to increase performance. Starting pipeline without parallelization!"

## HANDLE ARGUMENTS TODO: handle arguments better

assembly=$1
mnd=$2

pipeline=`cd "$( dirname $0)" && cd .. && pwd`
edit_mnd_script=${pipeline}/edit/edit-mnd-according-to-new-cprops.awk

[ $parallel == "true" ] && parallel -a ${mnd} --pipepart --will-cite --jobs 80% --block 1G "awk -f ${pipeline}/edit/edit-mnd-according-to-new-assembly.awk ${assembly} - " || awk -f ${pipeline}/edit/edit-mnd-according-to-new-assembly.awk ${assembly} ${mnd}
